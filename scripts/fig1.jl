using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor
using DataFrames, StatsPlots, DataFramesMeta, LaTeXStrings
using Distributions
### logistic

P = Dict{Symbol, Any}(
    :scaled => false,
    :S => 2:5:200,
    :μ => 1e-2,
    :σ => 5e-4,
    :k => 1.,
    :n0 => .01,
    :K => K,
    :dist => "normal",
    :N => 10,
    :seed => 10,
    :λ => 0.,
    :dist_r => Uniform(1.,2.)
);


r = []
r_se = []
for p in expand(P)
        stability!(p)
        push!(r, p[:richness])
        #push!(richness_se, p[:richness_se])
end

plot(
    P[:S], r, 
    #yerror = richness_se, 
    color = COLOR_LOG, 
    dpi = 500, 
    marker = :auto,
    label = "logistic",
    xlabel = "total number of species " * L"S",
    ylabel = "number of surviving species"
    )

### sublinear
P[:k] = .75
P[:K] = 1e6
r = []
r_se = []
for p in expand(P)
        stability!(p)
        push!(r, p[:richness])
        #push!(richness_se, p[:richness_se])
end

plot(
    P[:S], r, 
    #yerror = richness_se, 
    color = COLOR_SUB, 
    dpi = 500, 
    marker = :auto,
    label = "sublinear"
    )

savefig(plotsdir("surviving-species"))


##########

p = Dict{Symbol, Any}(
    :scaled => false,
    :S => 10,
    :μ => 1e-2,
    :λ => 0.,
    :σ => 5e-4,
    :k => 1.,
    :n0 => 1.,
    :K => K,
    :dist => "normal",
    :N => 5,
    :seed => rand(UInt),
    :dist_r => Uniform(.1,10)
);

evolve!(p, trajectory = true);
plot(p[:trajectory])