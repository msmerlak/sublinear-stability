using DrWatson
@quickactivate

using UnPack
using QuadGK, FastGaussQuadrature
using NLsolve
using LinearAlgebra

function H_RS(x, f, p)
    @unpack  ρ, μ, σ, β, V, λ = p
    n, z = x
    qd, q0, h = f
    @assert qd >= q0 "qd < q0"
    return - ρ^2*σ^2*β*(qd-q0)*n^2/2 + (ρ*μ*h-z*ρ*sqrt(q0)*σ)*n + V(n) #+ (1/β-λ)*log(n)
end

function F(f, p; n_min = 1e-2, n_max = 100, quadrature = gausshermite(3))
    @unpack  ρ, μ, σ, β, V, λ = p
    x, w = quadrature

    m(n, z) = exp(-β*H_RS([n, z], f, p))

    Z(z) = quadgk(n -> m(n, z), n_min, n_max)[1]
    N(z) = quadgk(n -> n*m(n, z), n_min, n_max)[1]
    N2(z) = quadgk(n -> n^2*m(n, z), n_min, n_max)[1]

    I1(z) = (1/sqrt(2π))*N2(z)/Z(z)
    I2(z) = (1/sqrt(2π))*(N(z)/Z(z))^2
    I3(z) = (1/sqrt(2π))*N(z)/Z(z)

    # qd, q0, h = F
    return [dot(w, Ik.(x)) for Ik in (I1, I2, I3)]
end


function F!(F, f, p; n_min = 1e-2, n_max = Inf, quadrature = gausshermite(3))
    @unpack  ρ, μ, σ, β, V, λ = p
    x, w = quadrature

    m(n, z) = exp(-β*H_RS([n, z], f, p))

    Z(z) = quadgk(n -> m(n, z), n_min, n_max)[1]
    N(z) = quadgk(n -> n*m(n, z), n_min, n_max)[1]
    N2(z) = quadgk(n -> n^2*m(n, z), n_min, n_max)[1]

    I1(z) = (1/sqrt(2π))*N2(z)/Z(z)
    I2(z) = (1/sqrt(2π))*(N(z)/Z(z))^2
    I3(z) = (1/sqrt(2π))*N(z)/Z(z)

    # qd, q0, h = F
    F[1] = dot(w, I1.(x))
    F[2] = dot(w, I2.(x))
    F[3] = dot(w, I3.(x))
end


function find_fp(p; n_min = 1e-3, n_max = 100, quadrature = gausshermite(3))
    @unpack  ρ, σ, β = p
    fp = fixedpoint((F, f) -> F!(F, f, p), 1e-2*[1., .5, 1.]; method = :anderson)
    qd, q0, h = fp.zero
    @assert qd >= q0 "ouch"
    return (qd, q0, h)
    #return (β*ρ*σ)^2(1 -  (β*ρ*σ)^2*(qd - q0)^2)
end


p = (
    ρ = 1.,
    β = 30.,
    μ = 2.,
    σ = 0.01,
    λ = 1.,
    V = n -> - (4/3)n^(3/4),
)

F(
    [0.41742, 0.35669, 0.58943],
    p;
    n_min = 0, n_max = 100, quadrature = gausshermite(5)    
)

find_fp(p)