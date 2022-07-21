using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))


using ProgressMeter, Suppressor, DataFrames
using Plots, LaTeXStrings

P = Dict{Symbol, Any}(
    :scaled => true,
    :S => 100,
    :μ => 1/100,
    :σ => 0.001,
    :k => 1.,
    :b0 => 1.,
    :K => 10,
    :λ => 0,
    :z => 0,
    :r => 1,
    #:dist_r => LogNormal(1,.1),
    :N => 1,
    :threshold => false,
    :dist => "normal",
    :symm => false,
    :seed => 736,
);

P[:K]=10^0
P[:r]=.54*P[:K]^(-1/4)

#P[:μ]*=P[:r]/P[:S] #mode(P[:dist_r])/P[:S]
#P[:σ]*=P[:μ]
#P[:K]*=P[:r]/(P[:σ]*sqrt(P[:S]+P[:μ]))
evolve!(P; trajectory = true);

#trajectories
plot(P[:trajectory],
ylabel = L"B_i(t)",
xlabel = L"t",
yscale = :log,
#yticks = [10^(-3),10^(-2),10^(-1), 10^(0), 10, 10^2, 10^3],
#ylims = [1e-4,1e0],
linewidth = 2,
legend = false,
alpha = .5,
grid = false,
palette = :blues, 
#palette = :YlOrBr_9,
)

#savefig("robustness-zero-mu5sigma01.svg")

# abunancies
histogram([(P[:equilibrium])],
normalize=true,
alpha = .5,
color = :gray,
ylabel = L"P(b^*)",
xlabel = L"b^* \textrm{(log10 \ scale)}",
label = false,
grid = false,
#ylims = [0,1.1],
#xlims = [-3,3],
#xticks = [-3,-2,-1,0,1,2,3],
)

plot!([i for i in 0:5], [P[:k]/(1+P[:K]*P[:S]*P[:μ]/P[:r]) for i in 0:5],
linewidth = 2,
linestyle = :dash,
linecolor = :gray,
)

X=[n for n in 0.1:.01:11]
#plot!((X),[pdf(Normal(P_n_log(P).μ, P_n_log(P).σ), n) for n in X],
plot!((X),[pdf(P_n_log(P), n) for n in X],
labels="cavity",
linewidth = 2,
linecolor = :black,
grid = false,
#xlims = [10^(-3),10^3],
#xticks = [10^-3,10^-2,10^-1,10^0,10^1,10^2,10^3],
)

# cavity solution with gaussian approximation
X=[n for n in .5*minimum(P[:equilibrium]):.0001:1.5*maximum(P[:equilibrium])]
plot!(X,[pdf(last(P_n_gauss(P)) ,n) for n in X],
labels="cavity - gaussian",
linewidth = 2,
alpha = 1,
linecolor = :red,
)

plot!(X,[P_n_mix(n, P) for n in X],
labels="cavity - mixed",
linewidth = 2,
alpha = 1,
linecolor = :black,
linestyle = :dash,
)

# cavity solution
ϕ, uno, due = Cavity(P, n_max=10^3)
X=[n for n in .5*minimum(P[:equilibrium]):.0001:1.5*maximum(P[:equilibrium])]
plot!(X,[P_n(n,uno,due,P) for n in X],
labels="cavity",
alpha = 1,
linewidth = 2,
linecolor = :black,
)

savefig("cavity+mixed+gaussian-s100.svg")

function f!(F, x)
    F[1] = 1-P[:S]*P_n(x[1], uno, due, P)/ϕ
end
n_max_new = nlsolve(f!, [1.]).zero[1]

# richness
P[:richness]/P[:S]
gauss_approx_ϕ(P, n_max=1000)

# spectrum
boundary(P, overprint=false)

# critical line
μ_critical, σ_critical = critical_line(P, μ_range=(.01:.4:1.21))
plot!(μ_critical, σ_critical,
labels = false,
linewidth = 4,
linecolor = COLOR_SUB35,
#ylims=[1.,20.],
#xlims=[.02,1.6]
)

plot!(μ_critical, .05*σ_critical.^0,
labels = false,
linewidth = 4,
linecolor = COLOR_LOG35,
#ylims=[0.01,.12],
#xlims=[.1,1.2]
)

# critical line for gaussian approximation
μ_critical, σ_critical = critical_line_mix(P, μ_range=(.1:.1:1.))
plot!(μ_critical, σ_critical,
labels = false,
linewidth = 4,
linecolor = COLOR_SUB35,
linestyle = :dash
#ylims=[1.,20.],
#xlims=[.02,1.6]
)

# critical line for gaussian approximation with alternative method
μ_critical, σ_critical = alternative_crit_σ(P, μ_range=(.01:.01:1.2))
plot!(μ_critical, σ_critical,
labels = false,
linewidth = 4,
linecolor = COLOR_SUB35,
#ylims=[1.,20.],
#xlims=[.02,1.6]
)

#savefig("papers/onofrio/fig/no-threshold-phase-space.svg")