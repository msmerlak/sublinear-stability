using DifferentialEquations
using Random, Distributions
using LinearAlgebra

function F(x, p)
    f = .- x.^2/p[:K]^(2 - p[:k]) .- x.*(p[:a]*x)

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
function evolve!(p)

    pb = ODEProblem(
        ODEFunction(
            (f, x, p, t) -> F!(f, x, p);
            jac = (j, x, p, t) -> J!(j, x, p)
            ),
            fill(5., p[:S]),
            (0., MAX_TIME),
            p
        )

    sol = solve(pb, 
        callback = CallbackSet(TerminateSteadyState(1e-3)), 
        save_on = false
        )
    p[:converged] = sol.retcode
    p[:richness] = mean(sol.u[end] .> 1e-3)
    p[:diversity] = Ω(sol.u[end])
end

function random_interactions!(p)

    (m, s) = p[:scaled] ? (p[:μ]/p[:S], p[:σ]/sqrt(p[:S])) : (p[:μ], p[:σ])

    if p[:dist] == "normal"
        dist = Normal(m - s, m + s)
    elseif p[:dist] == "uniform"
        dist = Uniform(min(0., m - s), max(2m, m + s))
    end

    a = rand(p[:rng], dist, (p[:S], p[:S]))
    a[diagind(a)] .= 0.

    p[:a] = a
end


function stats!(p)
    p[:rng] = MersenneTwister(p[:seed])

    stability = Vector{Bool}(undef, p[:N])
    richness = Vector{Float64}(undef, p[:N])
    diversity = Vector{Float64}(undef, p[:N])

    for i in 1:p[:N]
        random_interactions!(p)
        evolve!(p)
        stability[i] = p[:converged] == :Terminated
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