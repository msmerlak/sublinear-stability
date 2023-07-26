###

using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, ThreadsX
using Plots, LaTeXStrings, DelimitedFiles, Colors, ColorSchemes

# sublinear

P = Dict{Symbol, Any}(
        :scaled => false,
        :S => 2:2:100,
        :μ => .01,
        :σ => .005,
        :k => .75,
        :n0 => 1e-8,
        :λ => 0,
        :K => 1e6,
        :dist => "normal",
        #:dist_r => Uniform(.01,1),
        :N => 50,
        :symm => false,
    );

ϕ = ThreadsX.collect(full_coexistence(p) for p in expand(P));

open("../papers/onofrio/fig/SUB-stab-vs-div-N_$(P[:N])-S_$(P[:S])-n₀_$(P[:n0])-σ_$(P[:σ])-μ_$(P[:μ])-k_$(P[:k])-K_$(P[:K]).txt", "w") do io
    writedlm(io, [P[:S] ϕ], ',')
end 



# logistic

P = Dict{Symbol, Any}(
        :scaled => false,
        :S => 2:2:100,
        :μ => .01,
        :σ => .005,
        :k => 1.,
        :n0 => 1e-8,
        :λ => 0,
        :K => 20,
        :dist => "normal",
        #:dist_r => Uniform(.01,1),
        :N => 1e4,
        :symm => false,
    );

ϕ = ThreadsX.collect(full_coexistence(p) for p in expand(P));

open("../papers/onofrio/fig/LOG-stab-vs-div-N_$(P[:N])-S_$(P[:S])-n₀_$(P[:n0])-σ_$(P[:σ])-μ_$(P[:μ])-k_$(P[:k])-K_$(P[:K]).txt", "w") do io
    writedlm(io, [P[:S] ϕ], ',')
end


#=
a = readdlm("papers/onofrio/fig/stab-vs-div-N_100-R_100-S_2:2:100-n₀_1.0e-8-σ_0.005-μ_0.01-k_0.75-K_1.0e6.txt", ',')
b = readdlm("papers/onofrio/fig/stab-vs-div-N_100-R_100-S_2:2:100-n₀_1.0e-8-σ_0.005-μ_0.01-k_1.0-K_20.txt", ',')

plot(P[:S], b[:,1], 
ribbon=b[:,2],
c=COLOR_LOG49,
linewidth = 2,
grid = false,
label = "Logistic",
legend = :right,
)
plot!(P[:S], a[:,1], 
ribbon=a[:,2],
c=COLOR_SUB49,
linewidth = 2,
grid = false,
label = "Sublinear",
legend = :right,
)

plot!(P[:S], mean_ϕ, ribbon=std_ϕ,
c=COLOR_LOG49,
linewidth = 2,
legend = false
)
=#