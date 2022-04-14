using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, DataFrames
using Plots

P = Dict{Symbol, Any}(
    :scaled => false,
    :S => 50,
    :μ => .001,
    :σ => .00005,
    :k => .75,
    :n0 => 1,
    :K => 1e6,
    :dist => "normal",
    :N => 1,
    :seed => 1,
    :symm => false,
    :λ => 0,
    #:dist_r => LogNormal(1,.05),
);

evolve!(P; trajectory = true);
boundary(P, overprint=false)

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

uno, due = DMFT_eq(P, n_max=1000)
X=[n for n in minimum(P[:equilibrium]):.001:maximum(P[:equilibrium])]
plot!(X,[P_n(n,uno,due,P) for n in X],
labels=false,
linewidth = 2,
linecolor = :black,
)

# richness
P[:richness]/P[:S]
analytical_ϕ(P, uno, due, n_max=100)

# spectrum
boundary(P, overprint=false)

savefig("papers/onofrio/fig/logistic-spectra-Sfair.svg")

contour(X, Y, (x,y) -> T(x, y, n, p), levels = [1/(p[:σ])^2],
linewidth = 2, color = :auto, colorbar = false)