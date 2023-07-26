#= Produces the plots for the production scaling and for the macroecological laws. =#

using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, ThreadsX, DataFrames
using Plots, LaTeXStrings, DelimitedFiles, Colors, ColorSchemes, LinearAlgebra

#= Production statistics =#
P = Dict{Symbol, Any}(
        :scaled => false,
        :S => 50,
        :μ => 1, #10 .^(range(-2,stop=2,length=100)),
        :σ => 1.,
        :k => 1.,
        :b0 => 1.,
        :λ => 0,
        :K => 10 .^(range(-2,stop=2,length=5)),
        :dist => "normal",
        #:dist_r => Normal(.001,.0),
        :r => 1,
        :z => 0,
        :threshold => false,
        :N => 1,
        :symm => false,
        :seed => 34,

);

l = length(P[:K])

eq_distr = Matrix{Float64}(undef, P[:S], l)
growth_r = Matrix{Float64}(undef, P[:S], l)
carrying_cap = Matrix{Float64}(undef, P[:S], l)
cavity_moments = Matrix{Float64}(undef, 2, l)
richness = Vector{Float64}(undef, l)
for (i,p) in enumerate(expand(P))
    p[:r]=437*p[:K]^(-.25)
    p[:μ]=p[:r]/p[:K]/p[:S]
    #p[:μ]*=p[:r]/p[:S]/p[:b0]
    p[:σ]*=p[:μ]
    evolve!(p)
    richness[i] = p[:richness] 
    eq_distr[:, i] .= p[:equilibrium]
    growth_r[:, i] .= p[:r]
    carrying_cap[:, i] .= p[:K]
#    cavity_moments[1, i] = P_n_gauss(p)[1]
#    cavity_moments[2, i] = P_n_gauss(p)[2]
end

scatter([[mean(eq_distr[:,i]) for i in 1:l]],
[[mean(growth_r[:,i])*P[:b0]^(1-P[:k])*mean(eq_distr[:,i].^P[:k]) for i in 1:l]],
ylabel = L"\langle b* \ g(b^*)\rangle",
xlabel = L"\langle b*\rangle",
linewidth = 2,
label = false,
scale = :log,
alpha = .5,
xlims = [10^(-2),10^2],
#yticks = [10^(-2), 10^(-1),10^(0),10,10^(2)],
ylims=[1e-3,1e4]
)

open("production-sims.txt", "a") do io
        writedlm(io, [[mean(eq_distr[:,i]) for i in 1:l]  [mean(growth_r[:,i])*P[:b0]^(1-P[:k])*mean(eq_distr[:,i].^P[:k]) for i in 1:l]])
    end

vline!(1e-3:1e4, [P[:b0]*(P[:z]/P[:r])^(1/(P[:k]-1))],
color = :red, 
alpha = .5,
linewidth=2,
label = false,
linestyle = :solid)

# cavity - gaussian approximation
plot!([[cavity_moments[1,i] for i in 1:l]],
[mean(growth_r[:,i])*P[:b0]^(1-P[:k])*cavity_moments[1,i]^P[:k] for i in 1:l],
linewidth = 2,
linestyle = :dash,
alpha = .5,
linecolor = :black,
legends = :bottomright,
label = false ,
)

plot!([mean(eq_distr[:,i]) for i in 1:l],[(437/2^(2-.75))*mean(eq_distr[:,i])^.75 for i in 1:l],
linewidth = 2,
linestyle = :dash,
linecolor = :black,
alpha = .5,
legends = :bottomright,
#label = L"\langle b*\rangle^{3/4}" ,
label = false,
xlims = [10^(-2.5),10^2.5],
)   

#logistic
scatter!([[mean(eq_distr[:,i]) for i in 1:l]],
[[mean(growth_r[:,i])*mean(eq_distr[:,i] .* (1 .- eq_distr[:,i] ./ carrying_cap[:,i])) for i in 1:l]],
ylabel = L"\langle b* \ g(b^*)\rangle",
xlabel = L"\langle b*\rangle",
linewidth = 2,
label = false,
scale = :log,
alpha = .5,
xlims = [10^(-2.5),10^2.5],
#yticks = [10^(-2), 10^(-1),10^(0),10,10^(2)],
#ylims=[1e-6,1e-2]
#markersize = .1*richness,
grid = false,
)

#= Species abundance distributions =#
μᵣ = log(7.3)
σᵣ = .38

P = Dict{Symbol, Any}(
        :scaled => false,
        :S => 500,
        :μ => 1,
        :σ => .7,
        :k => .75,
        :b0 => .01,
        :λ => 0,
        :K => 1e10,
        :z => .1,
        :threshold => false,
        :dist => "normal",
        :dist_r => LogNormal(μᵣ,σᵣ),
        :N => 1,
        :symm => false,
        :seed => 21,
);

P[:μ]*=mode(P[:dist_r])/P[:S]
P[:σ]*=P[:μ]
P[:z]*=mode(P[:dist_r])

evolve!(P)

# plot for paper
# logscale histograms
histogram([log.(P[:equilibrium])],
normalize=true,
alpha = .5,
color = 2,
ylabel = L"P(N^*)",
xlabel = L"\log N^*",
label = false,
grid = false,
fillcolor= :red,
#ylims = [0,.75],
#yticks = [0,.1,.2,.3,.4,.5,.6,.7],
#xlims = [-5,5],
#xticks = [-5,-4,-3,-2,-1,0,1,2,3,4,5],
)

X=[n for n in -3:.001:3]
plot!((X),[pdf(Normal(P_n_log(P).μ, (P_n_log(P).σ)), n) for n in X],
labels="cavity",
linewidth = 2,
linecolor = :black,
grid = false,
)

histogram([P[:equilibrium]],
normalize=true,
ylabel = L"P(x)",
xlabel = L"x",
label = false,
grid = false,
ylims = [0,1.5],
#xlims=[0,maximum(log.(P[:equilibrium]))]
)
X=[n for n in 10^(-3):.001:20]
plot!((X),[pdf(P_n_log(P), n) for n in X],
labels="cavity",
linewidth = 2,
linecolor = :black,
grid = false,
)

# analytical approx for small σ (it works really good also for not so small)
X=[n for n in .5*minimum(P[:equilibrium]):.001:maximum(P[:equilibrium])]
plot!(X,[pdf(P_n_log(P), n) for n in X],
labels="cavity",
linewidth = 2,
linecolor = :black,
grid = false,
)
plot!(bn,[pdf(LogNormal(1,.1), n) for n in bn],
xscale=:log,
labels="cavity",
linewidth = 2,
linecolor = :black,
grid = false,
)

# logscale histograms
br=10 .^(log10(.75*minimum(P[:r])):.5:log10(5*maximum(P[:r])))
bn=10 .^(log10(.75*minimum(P[:equilibrium])):.5:log10(5*maximum(P[:equilibrium])))
hr = fit(Histogram{Float64}, P[:r], br)
hn = fit(Histogram{Float64}, P[:equilibrium], bn)
normalize!(hr; mode =:pdf) ; normalize!(hn; mode=:pdf)
hn.weights ./= diff(log10.(hn.edges[1]))

plot(hn; 
alpha=0.5, 
label = false,
xscale=:log10, 
xlim=[minimum(bn),maximum(bn)], 
grid=false,
xlabel = L"\log \ x",
ylabel = L"P(x)"
)

plot!(hr; 
alpha=0.5, 
label=[L"r" L"b^*"],
xscale=:log10, 
xlim=[minimum(min(br,bn)),maximum(max(br,bn))], 
grid=false,
xlabel = L"\log \ x",
ylabel = L"P(x)",
xticks = [10^(-3), 10^(-2), 10^(-1),10^(0), 10, 10^(2), 10^(3)],
)

#= Taylor's law =#
P = Dict{Symbol, Any}(
        :scaled => false,
        :S => 50,
        :μ => 10 .^(range(-4.9,stop=-.4,length=100)),
        :σ => 1.,
        :k => .75,
        :b0 => .01,
        :λ => 0,
        :K => 1e10,
        :r => .55,
        :z => 0,#490/10,
        :threshold => true,
        :dist => "gamma",
        #:dist_r => LogNormal(log(.1*(.01^.25)),.1),
        :N => 1,
        :symm => false,
        :seed => 31,
);

l = length(P[:μ])

eq_distr = Matrix{Float64}(undef, P[:S], l)
cavity = Matrix{Float64}(undef, 2, l)
for (i,p) in enumerate(expand(P))
    p[:μ]*=p[:r]/p[:S]/p[:b0] #mean(p[:dist_r])/p[:S]
    p[:σ]*=p[:μ]
    evolve!(p)
    eq_distr[:, i] = p[:equilibrium]
    cavity[1,i] = P_n_gauss(p)[1]
    cavity[2,i] = P_n_gauss(p)[2]-P_n_gauss(p)[1]^2
end

scatter([mean(eq_distr[:,i]) for i in 1:l], [var(eq_distr[:,i]) for i in 1:l],
scale = :log,
alpha = .5,
xlabel = L"\langle b*\rangle",
ylabel = L"\textrm{var}(b*)",
linewidth = 2,
label = false,
grid = false,
#color = :green,
)

P[:r]=4.67
P[:S]=50


ρ=1/(1-(1/(P[:k]-1))^2*P[:r]^2*P[:b0]^(2*(1-P[:k]))/P[:S])

plot!([mean(eq_distr[:,i]) for i in 1:100],[(ρ-1)*mean(eq_distr[:,i])^(2) for i in 1:100],
scale = :log,
linewidth = 2,
linestyle = :solid,
linecolor = :gray,
label = L"(\langle b\rangle^*)^{2}",
#yticks = [10^(-2), 10^(-1),10^(0),10,10^(2)],
#ylims=[1e-6,1e-2]
#xlims = [10^(-2),10^2],
legends = :topleft)

open("taylor-sims.txt", "a") do io
        writedlm(io, [[mean(eq_distr[:,i]) for i in 1:l]  [var(eq_distr[:,i]) for i in 1:l]])
    end

# cavity
plot!([cavity[1,i] for i in 1:5],[cavity[2,i] for i in 1:5],
scale = :log,
linewidth = 2,
linestyle = :dash,
linecolor = :gray,
label = L"(\langle b\rangle^*)^2",
#yticks = [10^(-2), 10^(-1),10^(0),10,10^(2)],
#ylims=[1e-2,1e-2],
legends = :topleft)

#= Size-density scaling =#
P = Dict{Symbol, Any}(
        :scaled => false,
        :S => 50,
        :μ => 10 .^(range(-4.8,stop=-.4,length=3)),
        :σ => 1,
        :k => .75,
        :b0 => .01,
        :λ => 0,
        :K => 1e10,
        :r => .55,
        :z => .55/10,
        :threshold => true,
        :dist => "gamma",
        #:dist_r => LogNormal(log(1/.1),log(.1)),
        :N => 1,
        :symm => false,
        :seed => 1,

);

l = length(P[:μ])

eq_distr = Matrix{Float64}(undef, P[:S], l)
growth_r = Matrix{Float64}(undef, P[:S], l)
for (i,p) in enumerate(expand(P))
    p[:μ]*=p[:r]/p[:S]/p[:b0]
    p[:σ]*=p[:μ]
    evolve!(p)
    eq_distr[:, i] = p[:equilibrium]
    growth_r[:, i] .= p[:r]
end

scatter!([mean((growth_r[:, i]/3.3).^(-4),dims=1) for i in 1:l],
[mean(eq_distr[:,i],dims=1)/mean((growth_r[:, i]/3.3).^(-4),dims=1) for i in 1:l],
ylabel = L"\langle n\rangle^*",
xlabel = L"m",
linewidth = 2,
label = false,
scale = :log,
ylims = [10^(-12),10^15],
yticks = [10^(-10),10^(-5),1,10^(5),10^10],
xlims = [10^(-12), 10^5],
xticks = [10^(-10),10^(-5),1,10^(5)],
)
plot!([eq_distr[:,i] for i in 1:l],[eq_distr[:,i].^(-1) for i in 1:l],
linewidth = 2,
linecolor = :black,
legends = :bottomright,
label = false,
)

open("size-density-sims.txt", "a") do io
        writedlm(io, [[mean((growth_r[:, i]/3.3).^(-4),dims=1) for i in 1:l]  [mean(eq_distr[:,i],dims=1)/mean((growth_r[:, i]/3.3).^(-4),dims=1) for i in 1:l]])
    end

scatter([[(growth_r[:, i]/(.01^.25)).^-4 for i in 1:l]],[[eq_distr[:,i] for i in 1:l]],
ylabel = L"\langle n^*\rangle",
xlabel = L"\langle m\rangle",
linewidth = 2,
label = false,
scale = :log,
)
plot!([eq_distr[:,i] for i in 1:l],[eq_distr[:,i].^(-1) for i in 1:l],
linewidth = 2,
linecolor = :black,
legends = :bottomright,
#label = L"\langle n_i^*\rangle^\alpha" ,
)