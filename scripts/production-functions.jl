using DrWatson
@quickactivate

include(srcdir("constants.jl"))

using Plots
using LaTeXStrings

K=100
plot(dpi = 500)
plot!(
    x -> (x - x^2/K)/((1-1/K)x),
    xlims = (0, K),
    xlabel = "Population density " * L"n",
    ylabel = "Per-capita productivity " * L"g(n)/n",
    label = "logistic",
    ylims = (1e-2,1),
    linewidth = 2
)
plot!(
    x -> x^.75/x,
    xlims = (0, K),
    color = COLOR_SUB,
    label = "sublinear",
    linewidth = 2
)
vline!([1], label = false)

savefig(plotsdir("productivity.pdf"))