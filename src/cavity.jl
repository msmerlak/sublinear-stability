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
        r = haskey(p, :r) ? first(p[:r]) : 1
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

# Expression for abundance with gaussian approximation
function P_n_gauss(p)
    @unpack μ, σ, k, S, b0 = p
    μ = !p[:scaled] ? μ*S : μ
    σ = !p[:scaled] ? σ*sqrt(S) : σ 
    r = haskey(p, :r) ? first(p[:r]) : 1
    e1 = (μ/r)^(1/(k-2))
    e2 = e1^2/(1-(1/(k-1))^2*(μ*e1/r)^(2*(2-k)/(k-1))*σ^2)
    ϕ = e2>0 ? (1+erf((e1-b0*p[:threshold])/sqrt(2*(e2-e1^2))))/2 : 1
    return (e1, e2>0 ? e2 : 0, ϕ, e2>0 ? Normal(e1,sqrt(e2-e1^2)) : Normal(e1,0))
end

# Expression for abundance mixed, complete solution with gaussian moments
function P_n_mix(n, p)
    @unpack μ, σ, k, S = p
    μ = !p[:scaled] ? μ*S : μ
    σ = !p[:scaled] ? σ*sqrt(S) : σ 
    r = haskey(p, :r) ? first(p[:r]) : 1
    e1 = (μ/r)^(1/(k-2))
    e2 = e1^2/(1-(1/(k-1))^2*(μ*e1/r)^(2*(2-k)/(k-1))*σ^2)
    return (1-k)*n^(k-2)/(sqrt(2*π*σ^2*e2/r^2))*
    exp(-(n^(k-1)-μ*e1/r)^2/(2*σ^2*e2/r^2))
end

# Expression for abundance distribution for uniform interaction and lognormal r
function P_n_log(p)
    μᵣ = p[:dist_r].μ
    σᵣ = p[:dist_r].σ
    μ = p[:scaled] ? p[:μ] : p[:μ]*p[:S]
    k = p[:k]

    e1 = exp((2*(1-k)*(μᵣ-log(μ/p[:b0]^(1-k)))+σᵣ^2)/(2*(1-k)*(2-k)))

    return LogNormal((μᵣ-log(μ*e1/p[:b0]^(1-k)))/(1-k), σᵣ/(1-k))
end

# Cavity calculation for fraction of survivals and first two moments of the equilibrium distribution given μ and σ.
function Cavity(p;
    Iter = 2000, #number of iteration
    rela = .01, #relax parameter for fixed point
    tol = 1e-9, #requested tolerance for numerical integrator
    n_max = 1e4, #upper bound for integration
    n_min = p[:b0]*p[:threshold], #lower bound for integration
    e1_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]), #initial guess for e1_n
    e2_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]) #initial guess for e1_n
    )

    ϕ_n = ones(Iter+1) #fraction of survivals
    e1_n = e1_init*ones(Iter+1) #<n>
    e2_n = e2_init*ones(Iter+1) #<n²>

    for i in 1:Iter
        ϕ_new = first(quadgk(x -> P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))
        e1_new = first(quadgk(x -> x*P_n(x, e1_n[i], e2_n[i], p)/ϕ_n[i] , n_min, n_max, rtol=tol))
        e2_new = first(quadgk(x -> x*x*P_n(x, e1_n[i], e2_n[i], p)/ϕ_n[i] , n_min, n_max, rtol=tol))

        ϕ_n[i+1] = (1-rela)*ϕ_n[i] + rela*ϕ_new
        e1_n[i+1] = (1-rela)*e1_n[i] + rela*e1_new
        e2_n[i+1] = (1-rela)*e2_n[i] + rela*e2_new
    end
    ϕ_n_eq, e1_n_eq, e2_n_eq = ϕ_n[Iter], e1_n[Iter], e2_n[Iter]
    return (ϕ_n_eq, e1_n_eq, e2_n_eq)
end

# Cavity calculation, variant through sampling of DMFT_eq. 
function Cavity_samples(p;
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

# Critical σ. Returns σ_cirit(μ), e1(μ,σ_crit), e1(μ,σ_crit)
function σ_crit(p;
    Iter = 2000, #number of iteration
    rela = .01, #relax parameter for fixed point
    tol = 1e-9, #requested tolerance for numerical integrator
    e1_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]), #initial guess for e1_n
    e2_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]), #initial guess for e1_n
    σc_init = (p[:scaled] ? .1*p[:μ] : .1*p[:μ]/p[:S]), #initial guess for σ_c
    ϵ = 1e-5, #Hadamard small parametrs for partie finie
    n_max = 1e3, #upper bound for integration
    n_min = p[:b0]*p[:threshold], #lower bound for integration
    )

    if !p[:scaled]  #singularity in ahmadian formula
        n_s = haskey(p, :r) ? (p[:μ]/(1-p[:k])/p[:r])^(1/(p[:k]-2)) : (p[:μ]/(1-p[:k]))^(1/(p[:k]-2))
    else
        n_s = haskey(p, :r) ? (p[:μ]/p[:S]/(1-p[:k])/p[:r])^(1/(p[:k]-2)) : (p[:μ]/p[:S]/(1-p[:k]))^(1/(p[:k]-2))
    end 

    ϕ_n = ones(Iter+1) 
    e1_n = e1_init*ones(Iter+1)
    e2_n = e2_init*ones(Iter+1) 
    σc = σc_init*ones(Iter+1) 
    
    if !haskey(p, :dist_r)
        for i in 1:Iter
            p[:σ] = σc[i]
            ϕ_new = first(quadgk(x -> P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))
            e1_new = first(quadgk(x -> x*P_n(x, e1_n[i], e2_n[i], p)/ϕ_n[i] , n_min, n_max, rtol=tol))
            e2_new = first(quadgk(x -> x*x*P_n(x, e1_n[i], e2_n[i], p)/ϕ_n[i] , n_min, n_max, rtol=tol))
            σc_new = ((first(quadgk(x -> P_n(x, e1_n[i], e2_n[i], p)/ϕ_n[i]/((1-p[:k])*x^(p[:k]-2)-p[:μ]/p[:S])^2 ,
            n_min, n_s - ϵ, rtol=tol)))/p[:S])^(-1/2)
            #σc_new = ((first(quadgk(x -> P_n(x, e1_n[i], e2_n[i], p)/ϕ_n[i]/((1-p[:k])*x^(p[:k]-2)-p[:μ]/p[:S])^2 , n_min, n_s - ϵ, rtol=tol))+
            #first(quadgk(x -> P_n(x, e1_n[i], e2_n[i], p)/ϕ_n[i]/((1-p[:k])*x^(p[:k]-2)-p[:μ]/p[:S])^2 , n_s + ϵ, n_max, rtol=tol))-
            #2*P_n(n_s, e1_n[i], e2_n[i], p)/ϵ)/p[:S])^(-1/2)
            
            ϕ_n[i+1] = (1-rela)*ϕ_n[i] + rela*ϕ_new
            e1_n[i+1] = (1-rela)*e1_n[i] + rela*e1_new
            e2_n[i+1] = (1-rela)*e2_n[i] + rela*e2_new
            σc[i+1] = (1-rela)*σc[i] + rela*σc_new
        end
    else
        a, b = p[:dist_r].a, p[:dist_r].b
        n_max = (p[:μ]/p[:S]/(1-p[:k])/b)^(1/(p[:k]-2)) - 1e-5
        for i in 1:Iter
            p[:σ] = σc[i]
            ϕ_new = first(quadgk(x -> P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))
            e1_new = first(quadgk(x -> x*P_n(x, e1_n[i], e2_n[i], p)/ϕ_n[i] , n_min, n_max, rtol=tol))
            e2_new = first(quadgk(x -> x*x*P_n(x, e1_n[i], e2_n[i], p)/ϕ_n[i] , n_min, n_max, rtol=tol))
            
            #hadamard regularization
            σc_new = (first(hcubature(x -> P_n(x[1], e1_n[i], e2_n[i], p)/ϕ_n[i]/(b-a)/(x[2]*(1-p[:k])*x[1]^(p[:k]-2)-p[:μ]/p[:S])^2, [n_min,n_max], [a,b], rtol=tol))/p[:S])^(-1/2)
        
            ϕ_n[i+1] = (1-rela)*ϕ_n[i] + rela*ϕ_new
            e1_n[i+1] = (1-rela)*e1_n[i] + rela*e1_new
            e2_n[i+1] = (1-rela)*e2_n[i] + rela*e2_new
            σc[i+1] = (1-rela)*σc[i] + rela*σc_new
            @show i
        end
    end
    σc_eq, ϕ_n_eq, e1_n_eq, e2_n_eq = σc[Iter], ϕ_n[Iter], e1_n[Iter], e2_n[Iter]
    return (σc_eq, ϕ_n_eq, e1_n_eq, e2_n_eq)
end

# Critical σ, variant with gaussian approximation
function σ_crit_gauss(p;
    Iter = 2000, #number of iteration
    rela = .01, #relax parameter for fixed point
    tol = 1e-9, #requested tolerance for numerical integrator
    ϵ = 1e-2, #Hadamard small parametrs for partie finie
    σc_init = (p[:scaled] ? .1*p[:μ] : .1*p[:μ]/p[:S]), #initial guess for σ_c
    n_max = 1e3, #upper bound for integration
    n_min = p[:b0]*p[:threshold], #lower bound for integration
    )

    if !p[:scaled]  #singularity in ahmadian formula
        n_s = haskey(p, :r) ? (p[:μ]/(1-p[:k])/p[:r])^(1/(p[:k]-2)) : (p[:μ]/(1-p[:k]))^(1/(p[:k]-2))
    else
        n_s = haskey(p, :r) ? (p[:μ]/p[:S]/(1-p[:k])/p[:r])^(1/(p[:k]-2)) : (p[:μ]/p[:S]/(1-p[:k]))^(1/(p[:k]-2))
    end 

    σc = σc_init*ones(Iter+1)     
    for i in 1:Iter
        p[:σ] = σc[i]
        σc_new = ((first(quadgk(x -> pdf(P_n_gauss(p),x)/((1-p[:k])*x^(p[:k]-2)-p[:μ]/p[:S])^2 ,
         n_min, n_s - ϵ, rtol=tol)))/p[:S])^(-1/2)

        σc[i+1] = (1-rela)*σc[i] + rela*σc_new
    end
    σc_eq = σc[Iter]
    return σc_eq
end

# Critical σ, variant with mixed approximation
function σ_crit_mix(p;
    Iter = 2000, #number of iteration
    rela = .01, #relax parameter for fixed point
    tol = 1e-9, #requested tolerance for numerical integrator
    ϵ = 1e-7, #Hadamard small parametrs for partie finie
    σc_init = (p[:scaled] ? .1*p[:μ] : .1*p[:μ]/p[:S]), #initial guess for σ_c
    n_min = p[:b0]*p[:threshold], #lower bound for integration
    )

    if !p[:scaled]  #singularity in ahmadian formula
        n_s = haskey(p, :r) ? (p[:μ]/(1-p[:k])/p[:r])^(1/(p[:k]-2)) : (p[:μ]/(1-p[:k]))^(1/(p[:k]-2))
    else
        n_s = haskey(p, :r) ? (p[:μ]/p[:S]/(1-p[:k])/p[:r])^(1/(p[:k]-2)) : (p[:μ]/p[:S]/(1-p[:k]))^(1/(p[:k]-2))
    end 

    σc = σc_init*ones(Iter+1)     
    for i in 1:Iter
        p[:σ] = σc[i]
        σc_new = ((first(quadgk(x -> P_n_mix(x,p)/((1-p[:k])*x^(p[:k]-2)-p[:μ]/p[:S])^2 ,
         n_min, n_s - ϵ, rtol=tol)))/p[:S])^(-1/2)

        σc[i+1] = (1-rela)*σc[i] + rela*σc_new
    end
    σc_eq = σc[Iter]
    return σc_eq
end

# Regularized critical σ. Returns σ_cirit(μ), e1(μ,σ_crit), e1(μ,σ_crit)
function σ_crit_reg(p;
    Iter = 2000, #number of iteration
    rela = .01, #relax parameter for fixed point
    tol = 1e-9, #requested tolerance for numerical integrator
    e1_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]), #initial guess for e1_n
    e2_init = (p[:scaled] ? 1/p[:μ] : 1/p[:μ]/p[:S]), #initial guess for e1_n
    σc_init = (p[:scaled] ? .1*p[:μ] : .1*p[:μ]/p[:S]), #initial guess for σ_c
    n_max = 5e3, #upper bound for integration
    n_min = p[:b0]*p[:threshold], #lower bound for integration
    )

    if !p[:scaled]  #singularity in ahmadian formula
        n_s = haskey(p, :r) ? (p[:μ]/(1-p[:k])/p[:r])^(1/(p[:k]-2)) : (p[:μ]/(1-p[:k]))^(1/(p[:k]-2))
    else
        n_s = haskey(p, :r) ? (p[:μ]/p[:S]/(1-p[:k])/p[:r])^(1/(p[:k]-2)) : (p[:μ]/p[:S]/(1-p[:k]))^(1/(p[:k]-2))
    end 

    ϕ_n = ones(Iter+1) 
    e1_n = e1_init*ones(Iter+1)
    e2_n = e2_init*ones(Iter+1) 
    σc = σc_init*ones(Iter+1) 
    
    for i in 1:Iter
        p[:σ] = σc[i]
        ϕ_new = first(quadgk(x -> P_n(x, e1_n[i], e2_n[i], p) , n_min, n_max, rtol=tol))
        e1_new = first(quadgk(x -> x*P_n(x, e1_n[i], e2_n[i], p)/ϕ_n[i] , n_min, n_max, rtol=tol))
        e2_new = first(quadgk(x -> x*x*P_n(x, e1_n[i], e2_n[i], p)/ϕ_n[i] , n_min, n_max, rtol=tol))
        σc_new = ((first(quadgk(x -> (ψ_reg(n_s+x, e1_n[i], e2_n[i], p, n_s)+
                 ψ_reg(n_s-x, e1_n[i], e2_n[i], p, n_s)-2*ψ_reg(n_s, e1_n[i], e2_n[i], p, n_s))/ϕ_n[i]/x , n_min, n_s, rtol=tol))+
                 (first(quadgk(x -> (ψ_reg(x, e1_n[i], e2_n[i], p, n_s)-
                 2*ψ_reg(n_s, e1_n[i], e2_n[i], p, n_s))/ϕ_n[i]/(x-n_s)^2 , 2*n_s, n_max, rtol=tol))))/p[:S])^(-1/2)
        
        ϕ_n[i+1] = (1-rela)*ϕ_n[i] + rela*ϕ_new
        e1_n[i+1] = (1-rela)*e1_n[i] + rela*e1_new
        e2_n[i+1] = (1-rela)*e2_n[i] + rela*e2_new
        σc[i+1] = (1-rela)*σc[i] + rela*σc_new
    end
    σc_eq, ϕ_n_eq, e1_n_eq, e2_n_eq = σc[Iter], ϕ_n[Iter], e1_n[Iter], e2_n[Iter]
    return (σc_eq, ϕ_n_eq, e1_n_eq, e2_n_eq)
end

function ψ_reg(n, e1_n, e2_n, p, n_s)
    r = haskey(p, :r) ? p[:r] : 1
    @unpack μ, σ, k, S = p
    P=(1-k)*n^(k-2)/(sqrt(2*π*σ^2*e2_n/r^2))*
    exp(-(n^(k-1)-μ*e1_n/r)^2/(2*σ^2*e2_n/r^2))
    if n == n_s
        return P*n_s^(2*(3-k))/r^2/(1-k)^2/(k-2)/(k-8)
    else
        return P*(n-n_s)^2/(1-k)^2/r^2/(n^(k-2)-n_s^(k-2))^2
    end
end


# critical line in μ - σ space
function critical_line(p;
    μ_range = (.005:.005:.1),
    )
    σc_line = ones(length(μ_range)) 
    for (i,μ) in enumerate(μ_range)
        p[:μ] = μ
        σc_line[i] = first(σ_crit(p))
        @show σc_line
    end
    return μ_range, σc_line
end

# critical line in μ - σ space for gaussian approximation
function critical_line_gauss(p;
    μ_range = (.005:.005:.1),
    )
    σc_line = ones(length(μ_range)) 
    for (i,μ) in enumerate(μ_range)
        p[:μ] = μ
        σc_line[i] = σ_crit_gauss(p)
        @show σc_line
    end
    return μ_range, σc_line
end

# critical line in μ - σ space for mixed approximation
function critical_line_mix(p;
    μ_range = (.005:.005:.1),
    )
    σc_line = ones(length(μ_range)) 
    for (i,μ) in enumerate(μ_range)
        p[:μ] = μ
        σc_line[i] = σ_crit_mix(p)
        @show σc_line
    end
    return μ_range, σc_line
end

# critical σ different approach
function alternative_crit_σ(p;
    μ_range = (.005:.005:.1),
    )
    σc_line = ones(length(μ_range)) 
    for (i,μ) in enumerate(μ_range)
        p[:μ] = μ
        @unpack μ, k, S = p
        μ = !p[:scaled] ? μ*S : μ
        r = haskey(p, :r) ? first(p[:r]) : 1
        e1 = (μ/r)^(1/(k-2))
        function f!(F, x)
            F[1] = 1/S-pdf(Normal(e1,sqrt(e1^2/(1-(1/(k-1))^2*(μ*e1/r)^(2*(2-k)/(k-1))*x[1]^2)-e1^2)), 0)
        end
        σc_line[i] = first(nlsolve(f!, [μ/10]).zero)
    end
    return μ_range, σc_line
end