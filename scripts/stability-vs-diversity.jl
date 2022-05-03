using Pkg
Pkg.add("DrWatson")

Pkg.instantiate()

using DrWatson
@quickactivate

using Glob
foreach(include, glob("*.jl", srcdir()))


using ProgressMeter, Suppressor, ThreadsX
using Plots, LaTeXStrings, DelimitedFiles, Colors, ColorSchemes

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
        :N => 10,
        :symm => false,
    );

R=10
meta_ϕ = Matrix{Float64}(undef, R, length(P[:S]))
for r in 1:R
    ϕ = []
    for p in expand(P) 
        append!(ϕ, full_coexistence(p))
    end
    meta_ϕ[r, :] = ϕ
end

mean_ϕ = vec(mean(meta_ϕ, dims=1))
std_ϕ = vec(std(meta_ϕ, dims=1))

plot!(P[:S], mean_ϕ, ribbon=std_ϕ,
c=COLOR_LOG49,
linewidth = 2,
legend = false
)

open("../papers/onofrio/fig/stab-vs-div-N_$(P[:N])-R_$R-S_$(P[:S])-n₀_$(P[:n0])-σ_$(P[:σ])-μ_$(P[:μ])-k_$(P[:k])-K_$(P[:K]).txt", "w") do io
    writedlm(io, [mean_ϕ std_ϕ], ',')
end 