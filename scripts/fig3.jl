using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor
using DataFrames, StatsPlots, DataFramesMeta, LaTeXStrings

### logistic

P = Dict{Symbol, Any}(
    :scaled => false,
    :S => 100:200:1000,
    :μ => 1e-3,
    :σ => 2e-4,
    :k => 1.,
    :n0 => 1.,
    :K => 100,
    :dist => "normal",
    :N => 1,
    :seed => rand(UInt)
);

plot(dpi = 500)
for p in expand(P)
        evolve!(p)
        density!(p[:equilibrium], label = "S = $(p[:S])", 
        xlabel = "abundance", ylabel = "density", linewidth = 2)
end
densities = current()


plot(dpi = 500)
for p in expand(P)
        boundary(p)
end
vline!([0.], linewidth = 2, color = "red", legend = false)
spectra = current()

logistic = plot(densities, spectra, layout = (2,1))

### sublinear
P[:k] = .75
P[:K] = 1e6

plot(dpi = 500)
for p in expand(P)
        evolve!(p)
        density!(p[:equilibrium], label = "S = $(p[:S])", 
        xlabel = "abundance", ylabel = "density", linewidth = 2)
end
densities = current()


plot(dpi = 500)
for p in expand(P)
        boundary(p)
end
#vline!([0.], linewidth = 2, color = "red", legend = false)
spectra = current()

sublinear = plot(densities, spectra, layout = (2,1))


### combine

plot(logistic)
savefig(plotsdir("spectra-logistic"))
plot(sublinear)
savefig(plotsdir("spectra-sublinear"))