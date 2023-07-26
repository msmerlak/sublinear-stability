#= Produces the plot of stability vs diversity (S) for fixed μ and σ. =#

using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, ThreadsX
using Plots, LaTeXStrings, DelimitedFiles, Colors, ColorSchemes

#= Sublinear =#
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

open("SUB-stab-vs-div-N_$(P[:N])-S_$(P[:S])-n₀_$(P[:n0])-σ_$(P[:σ])-μ_$(P[:μ])-k_$(P[:k])-K_$(P[:K]).txt", "w") do io
    writedlm(io, [P[:S] ϕ], ',')
end 

#= Logistic =#
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

open("LOG-stab-vs-div-N_$(P[:N])-S_$(P[:S])-n₀_$(P[:n0])-σ_$(P[:σ])-μ_$(P[:μ])-k_$(P[:k])-K_$(P[:K]).txt", "w") do io
    writedlm(io, [P[:S] ϕ], ',')
end

#= Plot =#
#=
a = readdlm("SUB ... .txt", ',')
b = readdlm("LOG ... .txt", ',')

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