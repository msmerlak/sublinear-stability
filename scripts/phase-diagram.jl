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
        :scaled => true,
        :S => 100,
        :μ => .1:.1:1.,
        :σ => .01:.01:.1,
        :k => .75,
        :n0 => 1e-8,
        :λ => 0,
        :K => 1e6,
        :dist => "normal",
        #:dist_r => Uniform(.01,1),
        :N => 5,
        :symm => false,
    );

# full coexistence

ϕ = ThreadsX.collect(full_coexistence(p) for p in expand(P));

open("papers/onofrio/fig/prova.txt", "w") do io
    writedlm(io, ϕ)
end

#prob-stability-S_$(P[:S])-N_$(P[:N])-n₀_$(P[:n0])-σ_$(P[:σ])-μ_$(P[:μ])

# sublinear = heatmap(
#     P[:μ], 
#     P[:σ],
#     reshape(ϕ, length(P[:μ]), length(P[:σ]))',
#     clims = (0,1),
#     legend = :none,
#     dpi = 500,
#     alpha = 1.,
#     c = palette([:white, COLOR_SUB49], 100),
#     grid = false,
#     xlabel = L"\mu = S\,\textrm{mean}(A_{ij})",
#     ylabel = L"\sigma = \sqrt{S}\,\textrm{sd}(A_{ij})",
#     title = L"\textrm{Probability \,\, of \,\, stability}"
# )

# vline!(P[:σ], [1],
# color = :black, 
# linewidth=2,
# linestyle = :dash)

# plot!(
#     μ -> critical_line_approx(μ, P[:k], P[:K], P[:S]), 
#     ylims = (minimum(P[:σ]), maximum(P[:σ])), 
#     linewidth = 2,
#     color = "blue"
#     )


# # diversity

# ϕ = ThreadsX.collect(diversity(p) for p in expand(P));
# sublinear = heatmap(
#     P[:μ], 
#     P[:σ],
#     reshape(ϕ, length(P[:μ]), length(P[:σ]))',
#     clims = (0,1),
#     legend = :none,
#     dpi = 500,
#     alpha = 1.,
#     c = palette([:white, COLOR_SUB49], 100),
#     grid = false,
#     xlabel = L"\mu = N\,\textrm{mean}(A_{ij})",
#     ylabel = L"\sigma = \sqrt{N}\,\textrm{sd}(A_{ij})",
#     title = "equilibrium diversity, S = $(P[:S]), K = $(P[:K]), $(P[:symm] ? "symmetric" : "non-symmetric")"
# )


# a = ThreadsX.collect(ahmadian(p) for p in expand(P));
# aa = heatmap(
#     P[:μ], 
#     P[:σ],
#     reshape(a, length(P[:μ]), length(P[:σ]))',
#     clims = (0,1),
#     legend = :none,
#     dpi = 500,
#     xlabel = L"\mu = N\,\textrm{mean}(A_{ij})",
#     ylabel = L"\sigma = \sqrt{N}\,\textrm{sd}(A_{ij})",
#     title = "ahmadian"
# )

# plot(sublinear, aa)

# ## scaled

# function critical_line_approx(μ, k, K, S)
#     if k == 1.
#         return (1/K - μ)
#     elseif k < 1.
#         return μ*(1-k-(2-k)/S)
#     end
# end

# ## unscaled

# function critical_line_approx(μ, k, K, S)
#     if k == 1.
#         return (1/K - μ)/sqrt(S)
#     elseif k < 1.
#         return μ*(1-k-(2-k)/S)*sqrt(S)
#     end
# end


# logistic = []
# sublinear = []
# for S in [10, 20, 30]

#     P = Dict{Symbol, Any}(
#         :scaled => false,
#         :S => S,
#         :μ => 1e-4:2e-3:2e-2,
#         :σ => 1e-4:2e-3:2e-2,
#         :k => 1.,
#         :n0 => 1e-5,
#         :λ => 0.,
#         :K => 100,
#         :dist => "normal",
#         #:dist_r => Uniform(1,2),
#         :N => 10,
#         :symm => true,
#         :seed => rand(UInt)
#     );


#     ϕ = ThreadsX.collect(diversity(p) for p in expand(P));
#     p = heatmap(
#         P[:μ], 
#         P[:σ],
#         reshape(ϕ, length(P[:μ]), length(P[:σ]))',
#         clims = (0,1),
#         legend = :none,
#         size = (800, 300)
#     )
#     plot!(
#         μ -> critical_line_approx(μ, P[:k], P[:K], P[:S]), 
#         ylims = (minimum(P[:σ]), maximum(P[:σ])), 
#         linewidth = 2,
#         color = "blue"
#         )
#     push!(logistic, p)

#     P[:k] = .75
#     P[:K] = 1e6
#     ϕ = ThreadsX.collect(diversity(p) for p in expand(P));

#     p = heatmap(
#         P[:μ], 
#         P[:σ],
#         reshape(ϕ, length(P[:μ]), length(P[:σ]))',
#         clims = (0,1),
#         legend = :none,
#         size = (800, 300)
#     )
#     plot!(
#         μ -> critical_line_approx(μ, P[:k], P[:K], P[:S]), 
#         ylims = (minimum(P[:σ]), maximum(P[:σ])), 
#         linewidth = 2,
#         color = "blue"
#         )
#     push!(sublinear, p)
# end

# unscaled = plot(logistic..., sublinear..., 
#     dpi = 500, size = (1000, 500)
#     )



# ### symmetric


# P = Dict{Symbol, Any}(
#     :scaled => true,
#     :S => 200,
#     :μ => .01:.2:9,
#     :σ => 0.01:.05:.6,
#     :k => .75,
#     :n0 => 1e-5,
#     :λ => 0.,
#     :K => 1e6,
#     :dist => "normal",
#     #:dist_r => Uniform(1,2),
#     :N => 50,
#     :symm => true
# );

# ϕ = ThreadsX.collect(diversity(p) for p in expand(P));
# sublinear_symmetric = heatmap(
# P[:μ], 
# P[:σ],
# reshape(ϕ, length(P[:μ]), length(P[:σ]))',
# clims = (0,1),
# legend = :none,
# dpi = 500,
# xlabel = L"\mu = N\,\textrm{mean}(A_{ij})",
# ylabel = L"\sigma = \sqrt{N}\,\textrm{sd}(A_{ij})",
# title = "equilibrium diversity"
# )

# =#
