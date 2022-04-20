using StatsBase
using Random, Distributions, Statistics
using DifferentialEquations
using LinearAlgebra
using QuadGK
using HCubature
using NLsolve
using SpecialFunctions

# Expression for abundance distribution in sublinear case
function P_n(n, e1_n, e2_n, p)
    @unpack μ, σ, k, S = p
    if !p[:scaled] μ = μ*S end
    if !p[:scaled] σ = σ*sqrt(S) end
    if !haskey(p, :dist_r) # unique r
        r = haskey(p, :r) ? p[:r] : 1
        return (1-k)*n^(k-2)/(sqrt(2*π*σ^2*e2_n/r^2))*
        exp(-(n^(k-1)-μ*e1_n/r)^2/(2*σ^2*e2_n/r^2))
    else # uniformely distributed r
        a, b = p[:dist_r].a, p[:dist_r].b
        return (1-k)*n^(k-2)/(sqrt(2*π*e2_n*σ^2)*(b-a)*2*n^(2*k-2))*
        sqrt(2*π*e2_n*σ^2)*((exp(-(a*n^(k-1)-μ*e1_n)^2/(2*σ^2*e2_n))-
        exp(-(b*n^(k-1)-μ*e1_n)^2/(2*σ^2*e2_n)))*sqrt(2*π*e2_n*σ^2)+
        sqrt(pi)*e1_n*μ*(erf((b*n^(k-1)-μ*e1_n)/sqrt(2*σ^2*e2_n))-
        erf((a*n^(k-1)-μ*e1_n)/sqrt(2*σ^2*e2_n))))
    end
end

# DMFT calculation for first two moments of the equilibrium distribution given μ and σ.
function DMFT_eq(p;
    Iter = 2000, #number of iteration
    rela = .01, #relax parameter for fixed point
    tol = 1e-9, #requested tolerance for numerical integrator
    n_max = (p[:μ]/p[:S]/(1-p[:k]))^(1/(p[:k]-2)) - 1e-5, #upper bound for integration
    n_min = p[:n0], #lower bound for integration
    e1_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]), #initial guess for e1_n
    e2_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]) #initial guess for e1_n
    )

    e1_n = e1_init*ones(Iter+1) #<n>
    e2_n = e2_init*ones(Iter+1) #<n²>

    for i in 1:Iter
        e1_new = first(quadgk(x -> x*P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))
        e2_new = first(quadgk(x -> x*x*P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))

        e1_n[i+1] = (1-rela)*e1_n[i] + rela*e1_new
        e2_n[i+1] = (1-rela)*e2_n[i] + rela*e2_new
    end
    e1_n_eq, e2_n_eq = e1_n[Iter], e2_n[Iter]
    return (e1_n_eq, e2_n_eq)
end

# DMFT calculation, variant through sampling of DMFT_eq. 
function DMFT_eq_samples(p;
    Iter = 2000, #number of iteration
    rela = .1, #relax parameter for fixed point
    N_samples = 10000, #number of samples
    e1_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]), #initial guess for e1_n
    e2_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]) #initial guess for e1_n
    )
    @unpack μ, σ, k, S = p
    if !p[:scaled] μ = μ*S end
    if !p[:scaled] σ = σ*sqrt(S) end

    e1_n = e1_init*ones(Iter+1) #<n>
    e2_n = e2_init*ones(Iter+1) #<n²>

    n = ones(Iter, N_samples)
    for i in 1:Iter
        η = rand(Normal(0,sqrt(e2_n[i])), N_samples)
        n[i,:] .= (μ*e1_n[i] .+ σ*η).^(1/(k-1))
    
        e1_n[i+1] = (1-rela)*e1_n[i] + rela*mean(n[i,:])
        e2_n[i+1] = (1-rela)*e2_n[i] + rela*mean(n[i,:].^2)
    end
    e1_n_eq, e2_n_eq = e1_n[Iter], e2_n[Iter]
    return (e1_n_eq, e2_n_eq)
end

# Analytical normalized richness
function analytical_ϕ(p, e1_n, e2_n;
    n_max = 1e4, 
    n_min = p[:n0],
    tol = 1e-9,
    )
    return first(quadgk(x -> P_n(x, e1_n, e2_n, p) , n_min, n_max, rtol=tol))
end

# Critical σ. Returns σ_cirit(μ), e1(μ,σ_crit), e1(μ,σ_crit)
function σ_crit(p;
    Iter = 1000, #number of iteration
    rela = .01, #relax parameter for fixed point
    tol = 1e-9, #requested tolerance for numerical integrator
    n_max = haskey(p, :r) ? (p[:μ]/p[:S]/(1-p[:k])/p[:r])^(1/(p[:k]-2)) - 1e-5 : (p[:μ]/p[:S]/(1-p[:k]))^(1/(p[:k]-2)) - 1e-5, #μ and α dependent upper cutoff for Ahmadian formula 
    n_min = p[:n0], #lower bound for integration
    e1_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]), #initial guess for e1_n
    e2_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]), #initial guess for e1_n
    σc_init = (p[:scaled] ? .1/p[:μ] : .1/p[:μ]/p[:S]), #initial guess for σ_c
    )

    e1_n = e1_init*ones(Iter+1)
    e2_n = e2_init*ones(Iter+1) 
    σc = σc_init*ones(Iter+1) 
    
    if !haskey(p, :dist_r)
        for i in 1:Iter
            p[:σ] = σc[i]
            e1_new = first(quadgk(x -> x*P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))
            e2_new = first(quadgk(x -> x*x*P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))
            σc_new = (first(quadgk(x -> P_n(x, e1_n[i], e2_n[i], p)/((1-p[:k])*x^(p[:k]-2)-p[:μ]/p[:S])^2 , n_min, n_max, rtol=tol))/p[:S])^(-1/2)

            e1_n[i+1] = (1-rela)*e1_n[i] + rela*e1_new
            e2_n[i+1] = (1-rela)*e2_n[i] + rela*e2_new
            σc[i+1] = (1-rela)*σc[i] + rela*σc_new
        end
    else
        a, b = p[:dist_r].a, p[:dist_r].b
        n_max = (p[:μ]/p[:S]/(1-p[:k])/b)^(1/(p[:k]-2)) - 1e-5
        for i in 1:Iter
            p[:σ] = σc[i]
            e1_new = first(quadgk(x -> x*P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))
            e2_new = first(quadgk(x -> x*x*P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))
            σc_new = (first(hcubature(x -> P_n(x[1], e1_n[i], e2_n[i], p)/(b-a)/(x[2]*(1-p[:k])*x[1]^(p[:k]-2)-p[:μ]/p[:S])^2, [n_min,n_max], [a,b], rtol=tol))/p[:S])^(-1/2)
            
            e1_n[i+1] = (1-rela)*e1_n[i] + rela*e1_new
            e2_n[i+1] = (1-rela)*e2_n[i] + rela*e2_new
            σc[i+1] = (1-rela)*σc[i] + rela*σc_new
            @show i
        end
    end
    σc_eq, e1_n_eq, e2_n_eq = σc[Iter], e1_n[Iter], e2_n[Iter]
    return (σc_eq, e1_n_eq, e2_n_eq)
end

# critical line in μ - σ space
function critical_line(p;
    μ_range = (.005:.005:.1),
    )
    σc_line = ones(length(μ_range)) 
    for (i,μ) in enumerate(μ_range)
        p[:μ] = μ
        σc_line[i] = first(σ_crit(p))
        #@show σc_line
    end
    return μ_range, σc_line
end