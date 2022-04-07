using DrWatson
@quickactivate

include(srcdir("constants.jl"))
using UnPack
using Plots
using LaTeXStrings


K=100

p = Dict{Symbol, Any}(
        :k => .75,
        :n0 => 1.,
        :K => 1e6
    );


plot(dpi = 500)

p[:k] = 1.
p[:K] = 100
plot(
    x -> production(x, p)/x,
    xlims = (0, K),
    color = COLOR_LOG,
    label = "logistic",
    linewidth = 2
)

p[:k] = .75
p[:K] = 1e6
plot!(
    x -> g(x, p)/x,
    xlims = (0, K),
    color = COLOR_SUB,
    label = "sublinear",
    linewidth = 2
)
vline!([p[:n0]], label = false)

savefig(plotsdir("productivity.pdf"))