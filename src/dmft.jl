using StatsBase
using Random, Distributions, Statistics
using DifferentialEquations
using LinearAlgebra
using QuadGK
using NLsolve

# Expression for abundance distribution in sublinear case
function P_n(n, e1_n, e2_n, p)
    @unpack μ, σ, k, S = p
    if !p[:scaled] μ = μ*S end
    if !p[:scaled] σ = σ*sqrt(S) end
    return (1-k)*n^(k-2)/(sqrt(2*π*e2_n*σ^2))*exp(-(n^(k-1)-μ*e1_n)^2/(2*σ^2*e2_n))
end

# DMFT calculation for first two moments of the equilibrium distribution
function DMFT_eq(p;
    Iter = 2000, #number of iteration
    a = .01, #relax parameter for fixed point
    tol = 1e-9, #requested tolerance for numerical integrator
    n_max = 1e2, #upper bound for integration
    n_min = p[:n0] #lower bound for integration
    )

    e1_n = ones(Iter+1) #<n>
    e2_n = ones(Iter+1) #<n²>

    for i in 1:Iter
        e1_new = first(quadgk(x -> x*P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))
        e2_new = first(quadgk(x -> x*x*P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))

        e1_n[i+1] = (1-a)*e1_n[i] + a*e1_new
        e2_n[i+1] = (1-a)*e2_n[i] + a*e2_new
    end
    e1_n_eq, e2_n_eq = e1_n[Iter], e2_n[Iter]
    return (e1_n_eq, e2_n_eq)
end

# Analytical normalized richness
function analytical_ϕ(p, e1_n, e2_n;
    n_max = 1e2, 
    n_min = p[:n0]
    )
    first(quadgk(x -> P_n(x, e1_n, e2_n, p) , n_min, n_max, rtol=tol))
end