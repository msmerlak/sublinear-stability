using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, DataFrames


P = Dict{Symbol, Any}(
    :scaled => false,
    :S => 50,
    :μ => 0.:.002:0.03,
    :σ => 0.:.002:0.03,
    :k => .75,
    :n0 => 1.,
    :K => 1e6,
    :dist => "normal",
    :N => 2,
    :seed => 1
);


# Threads.@threads for p in expand(P)
#     if p[:μ] > p[:σ]
#         name = savename(p, "jld2")
#         if !isfile(name)
#             stats!(p)
#             wsave(datadir("fig2", name), p)
#         end
#     end
# end


Threads.@threads for p in expand(P)
    name = savename(p, "jld2")
    if !isfile(name)
        stats!(p)
        wsave(datadir("fig2", name), p)
    end
end


