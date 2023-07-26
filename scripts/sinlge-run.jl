#= This script allows to evolve a community specified by P. It is possible to plot results and compare with predictions from the cavity solution. =#

using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, DataFrames
using Plots, LaTeXStrings

#= Define the system =#
P = Dict{Symbol, Any}(
    :scaled => true,
    :S => 10,
    :μ => 1.,
    :C => 1.,
    :σ => .05,
    :k => .75,
    :b0 => 1.,
    :x0 => 1.75,
    :K => 1e5,
    :λ => 0,
    :z => 0,
    :r => 1, 
    :N => 1,
    :threshold => false,
    :dist => "gamma",
    :symm => false,
    :seed => 17,
);

#= Evolve =#
evolve!(P; trajectory = true);

#= Trajectories =#
plot(P[:trajectory],
ylabel = L"B_i(t)",
xlabel = L"t",
#yscale = :log,
linewidth = 3,
legend = false,
alpha = 1,
grid = false,
palette = :YlOrBr_9,
)

P[:richness]

#= Density distribution =#
histogram([(P[:equilibrium])],
normalize=true,
alpha = .5,
color = :gray,
ylabel = L"P(B^*)",
xlabel = L"B^*",
label = false,
grid = false,
)

#= cavity solution with gaussian approximation =#
X=[n for n in .5*minimum(P[:equilibrium]):.0001:1.5*maximum(P[:equilibrium])]
plot!(X,[pdf(last(P_n_gauss(P)) ,n) for n in X],
labels="cavity - gaussian",
linewidth = 2,
alpha = 1,
linecolor = :black,
linestyle = :dash
)

#= cavity solution =#
ϕ, e1, e2 = Cavity(P, n_max=1000)
X=[n for n in .5*minimum(P[:equilibrium]):.0001:1.5*maximum(P[:equilibrium])]
plot!(X,[P_n(n,e1,e2,P) for n in X],
labels="cavity",
alpha = 1,
linewidth = 2,
linecolor = :black,
)

#= cavity solution with mixed gaussian approximation =#
plot!(X,[P_n_mix(n, P) for n in X],
labels = "cavity - mixed",
linewidth = 2,
alpha = 1,
linecolor = :black,
linestyle = :dashdot,
)

#= eigenvalues spectrum =#
boundary(P)