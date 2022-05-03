using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, ThreadsX, DataFrames
using Plots, LaTeXStrings, DelimitedFiles, Colors, ColorSchemes

### Production statistics ###

P = Dict{Symbol, Any}(
        :scaled => true,
        :S => 100,
        :μ => 10 .^(range(.1,stop=10,length=10)),
        :σ => .01,
        :k => .75,
        :n0 => 1e-8,
        :λ => 0,
        :K => 1e6,
        :dist => "normal",
        :dist_r => Normal(.005,.0),
        #:r => .005,
        :N => 1,
        :symm => false,
        :seed => 1,

);

l = length(P[:μ])

eq_distr = Matrix{Float64}(undef, P[:S], l)
growth_r = Matrix{Float64}(undef, P[:S], l)
for (i,p) in enumerate(expand(P))
    #p[:σ]*=p[:μ]
    evolve!(p)
    eq_distr[:, i] = p[:equilibrium]
    growth_r[:, i] = p[:r]
end

growth_r

scatter!([[mean(eq_distr[:,i]) for i in 1:l]],[[mean(growth_r.*eq_distr[:,i].^P[:k]) for i in 1:l]],
ylabel = L"\langle P(n^*)\rangle",
xlabel = L"\langle n^*\rangle",
linewidth = 2,
label = false,
scale = :log,
#ylims=[1e-6,1e-2]
)
plot!([mean(eq_distr[:,i]) for i in 1:l],[mean(eq_distr[:,i]).^P[:k] for i in 1:l],
linewidth = 2,
linecolor = :black,
legends = :bottomright,
label = L"\langle n_i^*\rangle^\alpha" ,
)   

### Species aboundance distributions ###

P = Dict{Symbol, Any}(
        :scaled => true,
        :S => 50,
        :μ => 1.,
        :σ => .05,
        :k => 1.,
        :n0 => .1,
        :λ => 0,
        :K => 1e1,
        :dist => "normal",
        :dist_r => Uniform(.00001,1),
        :N => 1,
        :symm => false,
        :seed => 2,
);

evolve!(P)
histogram([(P[:equilibrium]),(P[:r])],
normalize=true,
ylabel = "normalized species count",
xlabel = "population size",
label = false,
#xlims=[0,maximum(log.(P[:equilibrium]))]
)

# se riesci aggiungi linea analitica per lognormal

### Taylor's law ###

P = Dict{Symbol, Any}(
        :scaled => true,
        :S => 100,
        :μ => .1:.1:10.,
        :σ => .005,
        :k => .75,
        :n0 => 1e-8,
        :λ => 0,
        :K => 1e6,
        :dist => "normal",
        #:dist_r => Uniform(.01,1),
        :N => 1,
        :symm => false,
        :seed => 1,
        :r => 10,
);

l = length(P[:μ])

eq_distr = Matrix{Float64}(undef, P[:S], l)
for (i,p) in enumerate(expand(P))
    p[:σ]*=p[:μ]
    evolve!(p)
    eq_distr[:, i] = p[:equilibrium]
end

histogram(eq_distr, normalize = true)

scatter([mean(eq_distr[:,i]) for i in 1:l], [std(eq_distr[:,i]) for i in 1:l],
scale = :log,
xlabel = "mean(n*)",
ylabel = "std(n*)",
linewidth = 2,
label = false)
plot!([mean(eq_distr[:,i]) for i in 1:10],[.01mean(eq_distr[:,i])^.9 for i in 1:10],
scale = :log,
linewidth = 2,
linestyle = :dash,
linecolor = :black,
label = "linear",
legends = :topleft)

### Size-density scaling ###

P = Dict{Symbol, Any}(
        :scaled => true,
        :S => 100,
        :μ => 10 .^(range(.1,stop=1,length=10)),
        :σ => 1e-8,
        :k => .75,
        :n0 => 1e-8,
        :λ => 0,
        :K => 1e6,
        :dist => "normal",
        :dist_r => LogNormal(.5,.01),
        #:r => 1.,
        :N => 1,
        :symm => false,
        :seed => 1,

);

l = length(P[:μ])

eq_distr = Matrix{Float64}(undef, P[:S], l)
growth_r = Matrix{Float64}(undef, P[:S], l)
for (i,p) in enumerate(expand(P))
    p[:μ]*=mean(P[:dist_r])^(-4)
    p[:σ]*=p[:μ]
    evolve!(p)
    eq_distr[:, i] = p[:equilibrium]
    growth_r[:, i] = p[:r]
end

scatter!([[growth_r[:, i] for i in 1:l]],[[eq_distr[:,i] for i in 1:l]],
ylabel = L"\langle n^*\rangle",
xlabel = L"\langle r\rangle",
linewidth = 2,
label = false,
scale = :log,
)
plot!([eq_distr[:,i] for i in 1:l],[eq_distr[:,i].^P[:k] for i in 1:l],
linewidth = 2,
linecolor = :black,
legends = :bottomright,
label = L"\langle n_i^*\rangle^\alpha" ,
)

