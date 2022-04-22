using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, DataFrames
using Plots

P = Dict{Symbol, Any}(
    :scaled => true,
    :S => 1000,
    :μ => 1.,
    :σ => .1,
    :k => .75,
    :n0 => 1e-8,
    :K => 1e6,
    :dist => "normal",
    :N => 1,
    :seed => 1,
    :symm => false,
    :λ => 0,
    :r => 1,
    #:dist_r => Uniform(.1,1.),
);

evolve!(P; trajectory = true);

#trajectories
plot(P[:trajectory],
ylabel = "nᵢ(t)",
linewidth = 2,
legend = false,
)

# abunancies
histogram!([P[:equilibrium]],
normalize=true,
ylabel = "normalized species count",
xlabel = "population size",
label = false,
)

uno, due = DMFT_eq(P)
X=[n for n in .5*minimum(P[:equilibrium]):.0001:maximum(P[:equilibrium])]
plot!(X,[P_n(n,uno,due,P) for n in X],
labels=false,
linewidth = 2,
linecolor = :black,
)

# richness
P[:richness]/P[:S]
analytical_ϕ(P, uno, due, n_max=1000)

# spectrum
boundary(P, overprint=false)

# critical line
μ_critical, σ_critical = critical_line(P, μ_range=(.1:.1:1.))
plot!(μ_critical, σ_critical,
labels = false,
linewidth = 4,
linecolor = COLOR_SUB35,
fill = (0, 0.2, COLOR_SUB49),
#ylims=[1.,20.],
#xlims=[.02,1.6]
)
plot!(μ_critical, .03*σ_critical.^0,
labels = false,
linewidth = 4,
linecolor = COLOR_LOG35,
fill = (0, 0.2, COLOR_LOG49),
ylims=[0.01,.12],
xlims=[.1,1.2]
)


#savefig("papers/onofrio/fig/no-threshold-phase-space.svg")