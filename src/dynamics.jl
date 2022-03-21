using DifferentialEquations
using Random, Distributions
using LinearAlgebra

function F(x, p)
    # the vector field
    f = .- x.^2/p[:K]^(2 - p[:k]) .- x.*(p[:a]*x)

    # no growth at low densities, below n0
    above_threshold = x .> p[:n0]
    f[above_threshold] .+= x[above_threshold].^p[:k]
    return f
end

function F!(f, x, p)
    f .= .- x.^2/p[:K]^(2 - p[:k]) .- x.*(p[:a]*x)

    above_threshold = x .> p[:n0]
    f[above_threshold] .+= x[above_threshold].^p[:k]
end

function J(x, p)
    # J = F'
    j = - x.*p[:a]
    j[diagind(j)] .= .- 2x./p[:K]^(2 - p[:k]) .- a*x

    above_threshold = x .> p[:n0]
    j[diagind(j)[above_threshold]] .+= p[:k].*x[above_threshold].^(p[:k]-1)

    return j
end

function J!(j, x, p)
    j .= - x.*p[:a]
    j[diagind(j)] .= .- 2x./p[:K]^(2 - p[:k]) .- p[:a]*x

    above_threshold = x .> p[:n0]
    j[diagind(j)[above_threshold]] .+= p[:k].*x[above_threshold].^(p[:k]-1)
end


## solving

MAX_TIME = 10_000

invasion() = DiscreteCallback((u, t, integrator) -> Ω(u) < .05, terminate!)

function evolve!(p)

    pb = ODEProblem(
        ODEFunction(
            (f, x, p, t) -> F!(f, x, p); #in-place F faster
            jac = (j, x, p, t) -> J!(j, x, p) #specify jacobian speeds things up
            ),
            fill(5., p[:S]), #initial condition
            (0., MAX_TIME),
            p
        )

    sol = solve(pb, 
        callback = CallbackSet(TerminateSteadyState(1e-3), invasion()), 
        save_on = false #don't save whole trajectory, only endpoint
        )
    p[:converged] = (sol.retcode == :Terminated && Ω(sol.u[end]) > .05)
    p[:richness] = mean(sol.u[end] .> p[:n0])
    p[:diversity] = Ω(sol.u[end])
end

function random_interactions!(p)
    # add a random interaction matrix to p, the dict of parameters
    (m, s) = p[:scaled] ? (p[:μ]/p[:S], p[:σ]/sqrt(p[:S])) : (p[:μ], p[:σ])

    if p[:dist] == "normal"
        dist = Normal(m, s)
    elseif p[:dist] == "uniform"
        dist = Uniform(min(0., m - s), max(2m, m + s))
    end

    a = rand(p[:rng], dist, (p[:S], p[:S]))
    a[diagind(a)] .= 0. #self-regulation is not part of interaction matrix

    p[:a] = a
end


function stats!(p)
    # run N simulates and append results to p
    p[:rng] = MersenneTwister(p[:seed])

    stability = Vector{Bool}(undef, p[:N])
    richness = Vector{Float64}(undef, p[:N])
    diversity = Vector{Float64}(undef, p[:N])

    for i in 1:p[:N]
        random_interactions!(p)
        evolve!(p)
        stability[i] = p[:converged]
        richness[i] = p[:richness]
        diversity[i] = p[:diversity]
    end
    delete!(p, :converged)
    delete!(p, :equilibrium)
    delete!(p, :rng)

    p[:prob_stab] = mean(stability)

    p[:richness] = mean(richness)
    p[:richness_se] = std(richness)/sqrt(p[:N])

    p[:diversity] = mean(diversity)
    p[:diversity_se] = std(diversity)/sqrt(p[:N])
end