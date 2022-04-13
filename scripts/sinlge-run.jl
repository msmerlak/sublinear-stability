using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, DataFrames
using Plots

P = Dict{Symbol, Any}(
    :scaled => false,
    :S => 50,
    :μ => .001,
    :σ => .0001,
    :k => 1.,
    :n0 => 1,
    :K => 100,
    :dist => "normal",
    :N => 1,
    :seed => 1,
    :symm => false,
    :λ => 0,
    #:dist_r => LogNormal(.5,.25),
);

evolve!(P; trajectory = true);

#trajectories
plot(P[:trajectory],
ylabel = "nᵢ(t)",
linewidth = 2,
legend = false,
)

# abunancies
histogram([P[:equilibrium]],
normalize=true,
ylabel = "normalized species count",
xlabel = "population size",
label = false,
)

P[:richness]/P[:S]

# spectrum
boundary(P)

savefig("logistic-spectra.svg")