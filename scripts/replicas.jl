using DrWatson
@quickactivate

using UnPack
using QuadGK, FastGaussQuadrature
using NLsolve, SpeedMapping
using LinearAlgebra

function H_RS(x, f, p)
    @unpack  μ, σ, β, V= p
    n, z = x
    qd, q0, h = f
    return - σ^2*β*(qd-q0)*n^2/2 + (μ*h-z*sqrt(q0)*σ)*n + V(n) #+ (1/β-λ)*log(n)
end

function int(f, a, b)
    x, w = gausslegendre(5)
    d = (b-a)/2
    s = (b+a)/2
    return d*dot(w, f.(d*x .+ s)) 
end

function F(f, p; n_min = 1e-2, n_max = 100, quadrature = gausshermite(3))

    @unpack  μ, σ, β, V= p
    x, w = quadrature

    m(n, z) = exp(-β*H_RS([n, z], f, p))

    # Z(z) = quadgk(n -> m(n, z), n_min, n_max)[1]
    # N(z) = quadgk(n -> n*m(n, z), n_min, n_max)[1]
    # N2(z) = quadgk(n -> n^2*m(n, z), n_min, n_max)[1]

    Z(z) = int(n -> m(n, z), n_min, n_max)
    N(z) = int(n -> n*m(n, z), n_min, n_max)
    N2(z) = int(n -> n^2*m(n, z), n_min, n_max)

    println(Z(.1))
    I1(z) = N2(z*sqrt(2))/Z(z*sqrt(2))
    I2(z) = (N(z*sqrt(2))/Z(z*sqrt(2)))^2
    I3(z) = N(z*sqrt(2))/Z(z*sqrt(2))

    # qd, q0, h = F
    return (1/sqrt(π))*[dot(w, Ik.(x)) for Ik in (I1, I2, I3)]
end


plot(n -> H_RS([n, .1], [1., .2, .1], p))
F([1., .2, .1], p; n_min = 1., n_max = 3)

function replicon(p; kwargs...)

    @unpack  σ, β = p


    fp = fixedpoint(f -> F(f, p; kwargs...), [.5, .4, .5]; method = :anderson, m = 10, beta = .1, store_trace = true)
    
    @assert fp.f_converged "FP iteration didn't converge"

    qd, q0, h = fp.zero
    @show [qd, q0, h]
    return (β*σ)^2(1 - (β*σ)^2*(qd - q0)^2)
end

p = (
    β = 50.,
    μ = 5.,
    σ = 0.3,
    V = n -> - (4/3)n^(3/4),
)

F([.5, .1, .1], p)

replicon(p; n_max = 100)

# Ada's numbers = [0.41742, 0.35669, 0.58943]
