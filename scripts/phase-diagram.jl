#= This script allows to produce plots for stability in the μS - σ√S plane. =#

using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, ThreadsX
using Plots, LaTeXStrings, DelimitedFiles, Colors, ColorSchemes

#= Define P =#
P = Dict{Symbol, Any}(
        :scaled => true,
        :S => 100,
        :μ => .2:.2:10.,
        :σ => .1:.1:5.,
        :k => .75,
        :b0 => 1.,
        :λ => 0,
        :z => 1.,
        :K => 2,
        :threshold => false,
        :dist => "normal",
        #:dist_r => Uniform(.01,1),
        :N => 5,
        :symm => false,
    );

#= Full stable coexistence =#
ϕ = ThreadsX.collect(full_coexistence(p) for p in expand(P));

#= Write down results =#
a = reshape(ϕ, length(P[:μ]), length(P[:σ]))
open("robustness-carrying-c-fixed-S_$(P[:S])-N_$(P[:N])-σ_$(P[:σ])-μ_$(P[:μ])-K_$(P[:K]).txt", "w") do io
    writedlm(io, a)
end

#= Plot heatmap =#
sublinear = heatmap(
P[:μ], 
P[:σ],
a',
clims = (0,1),
legend = :none,
dpi = 500,
alpha = 1.,
c = palette([:white, COLOR_SUB49], 100),
grid = false,
xlabel = L"\mu S",
ylabel = L"\sigma \sqrt{S}",
title = L"\textrm{Stability \,\, in \,\, parameter \,\, space}"
)

#= critical line =#
μ_critical, σ_critical = critical_line(P, μ_range=(.01:.1:1.21))
plot!(μ_critical, σ_critical,
labels = false,
linewidth = 4,
linecolor = COLOR_SUB35,
ylims=[.00,.117],
#xlims=[.02,1.]
)

plot!(μ_critical, .05*σ_critical.^0,
labels = false,
linewidth = 4,
linecolor = COLOR_LOG35,
#ylims=[0.01,.12],
#xlims=[.1,1.2]
)

#= critical line for gaussian approximation =#
μ_critical, σ_critical = critical_line_gauss(P, μ_range=(.01:.4:1.21))
plot!(μ_critical, σ_critical,
labels = false,
linewidth = 4,
linecolor = COLOR_LOG35,
linestyle = :dash,
ylims=[.0,.117],
#xlims=[.02,1.6]
)

#= Feasibility threshold =#
ν = sqrt(2*log(P[:S]))

n = [x for x in .01:.005:1]
plot(n, ((1-P[:k])/ν)*(n .- n.^((P[:k]-3)/(P[:k]-2))),
color = :green,
linestyle = :solid,
linewidth = 4,
alpha=1,
ylims = [0.001,0.04],
xlims = [.0,1]
)

n = [x for x in .01:.005:1.2]
plot!(n, ((1-P[:k])/ν)*n,
color = :green,
linestyle = :solid,
linewidth = 4,
alpha=1,
ylims = [0.,.12],
#xlims = [.0,1]
)