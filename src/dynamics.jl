using DifferentialEquations
using Random, Distributions
using LinearAlgebra
using Measurements

function production(x, p)
        return (x > p[:b0]*p[:threshold] ? x^p[:k] : 0) - x^2/p[:K]
end
function dproduction(x, p)
    return (x > p[:b0]*p[:threshold] ? p[:k]*x^(p[:k]-1) : 0) - 2x/p[:K]
end

function F!(f, x, p)
    f .= p[:r]*p[:b0]^(1-p[:k]).*production.(ppart.(x), Ref(p)) .- p[:z]*x .- x.*(p[:a]*x) .+ p[:λ]
end

function J!(j, x, p)
    j = - x.*p[:a]
    j[diagind(j)] .= p[:r]*p[:b0]^(1-p[:k]).*dproduction.(x, Ref(p)) .- 2x./p[:K]^(2 - p[:k]) .-p[:z] .- p[:a]*x


    # above_threshold = x .> p[:b0]
    # j[diagind(j)[above_threshold]] .+= p[:k].*p[:r][above_threshold].*x[above_threshold].^(p[:k]-1)
end


## solving

MAX_TIME = 1e3
MAX_ABUNDANCE = 1e5

blowup() = DiscreteCallback((u, t, integrator) -> maximum(u) > MAX_ABUNDANCE, terminate!)
converged() = TerminateSteadyState(1e-4)

function evolve!(p; trajectory = false)

    if !haskey(p, :rng) p[:rng] = MersenneTwister(p[:seed]) end
    if !haskey(p, :a) add_interactions!(p) end
    if !haskey(p, :r) add_growth_rates!(p) end
    if !haskey(p, :x0) add_initial_condition!(p) end

    pb = ODEProblem(
        ODEFunction(
            (f, x, p, t) -> F!(f, x, p); #in-place F faster
            jac = (j, x, p, t) -> J!(j, x, p) #specify jacobian speeds things up
            ),
            fill(1., p[:S]), #initial condition
            (0., MAX_TIME),
            p
        )

    sol = solve(pb, 
        callback = CallbackSet(converged(), blowup()), 
        save_on = trajectory #don't save whole trajectory, only endpoint
        )

    p[:equilibrium] = sol.retcode == :Terminated ? sol.u[end] : NaN
    p[:converged] = (sol.retcode == :Terminated && maximum(p[:equilibrium]) < MAX_ABUNDANCE)
    p[:richness] = sum(sol.u[end] .> p[:b0]*p[:threshold])
    p[:diversity] = p[:richness] == 0 ? 0 : Ω(sol.u[end].*(sol.u[end] .> p[:b0]*p[:threshold]))

    p[:survivors] = count(>=(p[:b0] * p[:threshold]), sol.u[end])
    p[:biomass] = sum(sol.u[end])

    if trajectory
        p[:trajectory] = sol
    end

end

function add_initial_condition!(p)
    p[:x0] = rand(p[:rng], Uniform(2, 10), p[:S])
end

function equilibria!(p)
    if !haskey(p, :rng) p[:rng] = MersenneTwister(p[:seed]) end
    if !haskey(p, :a) add_interactions!(p) end
    if !haskey(p, :r) add_growth_rates!(p) end

    equilibria = Vector{Float64}[]
    sizehint!(equilibria, p[:N])

    Threads.@threads for _ in 1:p[:N]
        add_initial_condition!(p)
        pb = ODEProblem(
            ODEFunction(
                (f, x, p, t) -> F!(f, x, p); #in-place F faster
                jac = (j, x, p, t) -> J!(j, x, p) #specify jacobian speeds things up
                ),
                p[:x0],
                (0., MAX_TIME),
                p
            )
        sol = solve(pb, 
            callback = CallbackSet(TerminateSteadyState(1e-5), blowup()), 
            save_on = false #don't save whole trajectory, only endpoint
            )
        push!(equilibria, sol.u[end])
    end
    p[:equilibria] = uniquetol(equilibria, atol = .1)
    p[:num_equilibria] = length(p[:equilibria])
    p[:num_interior_equilibria] = sum(map(x-> all(x .> p[:b0]*p[:threshold]), p[:equilibria]))

    return p[:num_equilibria]
end



function add_interactions!(p)
    # add a random interaction matrix to p, the dict of parameters
    (m, s) = p[:scaled] ? (p[:μ]/p[:S], p[:σ]/sqrt(p[:S])) : (p[:μ], p[:σ])

    if p[:dist] == "normal"
        dist = Normal(m, s)
    elseif p[:dist] == "uniform"
        dist = Uniform(max(0., m - s), min(2m, m + s))
    elseif p[:dist] == "gamma"
        dist = Gamma(m^2/s^2, s^2/m)
    end

    a = rand(p[:rng], dist, (p[:S], p[:S]))
    a[diagind(a)] .= 0. #self-regulation is not part of interaction matrix


    if p[:symm] 
        for i in 1:p[:S], j in 1:p[:S]
            if i > j
                a[i, j] = a[j, i]
            end
        end
    end
    p[:a] = a
end

function add_growth_rates!(p)
    if haskey(p, :dist_r)
        p[:r] = rand(p[:rng], p[:dist_r], p[:S])
    else 
        p[:r] = ones(p[:S])
    end
end

function stability!(p)
    # run N simulates and append results to p

    p[:rng] = MersenneTwister(p[:seed])
    stability = Vector{Bool}(undef, p[:N])
    diversity = Vector{Float64}(undef, p[:N])
    richness = Vector{Float64}(undef, p[:N])

    Threads.@threads for i in 1:p[:N]
        add_interactions!(p)
        add_growth_rates!(p)
        evolve!(p)
        stability[i] = p[:converged]
        diversity[i] = p[:diversity]
        richness[i] = p[:richness]
    end
    delete!(p, :converged)
    delete!(p, :rng)

    p[:richness] = mean(richness)/p[:S]
    p[:prob_stab] = mean(stability)


    p[:diversity] = mean(diversity)
    p[:diversity_se] = std(diversity)/sqrt(p[:N])
end




function diversity(p)
    # run N simulates and append results to p
    diversity = Vector{Float64}(undef, p[:N])
    p[:rng] = MersenneTwister()

    for i in 1:p[:N]
        add_interactions!(p)
        add_growth_rates!(p)
        add_initial_condition!(p)
        evolve!(p)
        diversity[i] = p[:diversity]
    end

    return mean(diversity)
end

function full_coexistence(p)
    # run N simulates and append results to p
    full_coexistence = Vector{Float64}(undef, p[:N])
    p[:rng] = MersenneTwister()

    for i in 1:p[:N]
        add_interactions!(p)
        add_growth_rates!(p)
        add_initial_condition!(p)
        evolve!(p)
        full_coexistence[i] = p[:converged] && p[:richness]==p[:S] ? 1 : 0
    end

    return mean(full_coexistence)
end

function survivors(p)
    # run N simulates and append results to p
    survivors = Vector{Float64}(undef, p[:N])
    p[:rng] = MersenneTwister()

    for i in 1:p[:N]
        add_interactions!(p)
        add_growth_rates!(p)
        add_initial_condition!(p)
        evolve!(p)
        survivors[i] = p[:survivors]
    end

    return mean(survivors) .± std(survivors)/sqrt(p[:N])
end

function biomass(p)
    # run N simulates and append results to p
    biomass = Vector{Float64}(undef, p[:N])
    p[:rng] = MersenneTwister()

    for i in 1:p[:N]
        add_interactions!(p)
        add_growth_rates!(p)
        add_initial_condition!(p)
        evolve!(p)
        biomass[i] = p[:biomass]
    end

    return mean(biomass) .± std(biomass)/sqrt(p[:N])
end

function ahmadian(p)
    # run N simulates and append results to p
    ahmadian = Vector{Bool}(undef, p[:N])
    p[:rng] = MersenneTwister()
    @unpack μ, σ, k = p

    for i in 1:p[:N]
        add_interactions!(p)
        add_growth_rates!(p)
        add_initial_condition!(p)
        evolve!(p)
        ahmadian[i] = sum(1 ./(μ .- (1-k)*ppart.(p[:equilibrium]).^(k-2))) < 1/σ^2
    end

    return mean(ahmadian)
end