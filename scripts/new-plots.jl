using DrWatson, Glob, Revise
@quickactivate
foreach(includet, glob("*.jl", srcdir()))

using Plots; gr(dpi = 500)
using LaTeXStrings

P(k, S, K, μ, r = 1, N = 100) = Dict{Symbol, Any}(
    :scaled => false,
    :threshold => 1,
    :b0 => 0.01,
    :r => 1,
    :z => .1,
    :S => S,
    :μ => μ,
    :σ => μ,
    :k => k,
    :λ => 0,
    :K => k == 1 ? K : Inf,
    :dist => "gamma",
    :N => N,
    :seed => rand(UInt),
    :symm => false
);


# vary S
S_range = 5:5:100
z = .1
for μ ∈ (.01, .1, 1.)

    plot(
        xlabel = "Number of species " * L"S",
        ylabel = "Fraction of surviving species",
        title = "r = 1, z = $z, μ = σ = $μ (gamma)"
    )

    k = 1
    for K ∈ (.1, 100) 
        plot!(
            S_range,
            survivors.(P.(k, S_range, K, μ)) ./ S_range,
            label = "k = $k, K = $K",
            markers = :auto,
            color = k == 1 ? COLOR_LOG : COLOR_SUB,
        )
    end 

    k = .75
    plot!(
        S_range,
        survivors.(P.(k, S_range, K, μ)) ./ S_range,
        label = "k = $k",
        color = k == 1 ? COLOR_LOG : COLOR_SUB,
    )
end
current()

p = P(.75, 50, Inf, 5., 1., 1)
evolve!(p; trajectory = true)
plot(
    p[:trajectory],
    legend = false
)