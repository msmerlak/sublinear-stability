using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor
using DataFrames, StatsPlots, DataFramesMeta, LaTeXStrings

### logistic

P = Dict{Symbol, Any}(
    :scaled => false,
    :S => 2:1:50,
    :μ => 1e-2,
    :σ => 2e-4,
    :k => 1.,
    :n0 => 1.,
    :K => K,
    :dist => "normal",
    :N => 100,
    :seed => rand(UInt)
);


richness = []
richness_se = []
for p in expand(P)
        stability!(p)
        push!(richness, p[:richness])
        push!(richness_se, p[:richness_se])
end

plot(
    P[:S], richness, 
    yerror = richness_se, 
    color = COLOR_LOG, 
    dpi = 500, 
    marker = :auto,
    label = "logistic",
    xlabel = "total number of species " * L"S",
    ylabel = "fraction of surviving species"
    )

### sublinear
P[:k] = .75
P[:K] = 1e6

richness = []
richness_se = []
for p in expand(P)
        stability!(p)
        push!(richness, p[:richness])
        push!(richness_se, p[:richness_se])
end

plot!(
    P[:S], richness, 
    yerror = richness_se, 
    color = COLOR_SUB, 
    dpi = 500, 
    marker = :auto,
    label = "sublinear"
    )

savefig(plotsdir("surviving-species"))