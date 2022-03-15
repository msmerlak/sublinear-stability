using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor

@suppress_err begin 
    P = Dict{Symbol, Any}(
        :scaled => false,
        :S => 5:5:150,
        :μ => 5e-3,
        :σ => 1e-4,
        :k => [.75, 1.],
        :n0 => [0., 1.],
        :K => [20, 50, 100, 500, 1e10],
        :dist => "uniform",
        :N => 10,
        :seed => 1
    );

    Threads.@threads for p in expand(P)
        name = savename(p, "jld2")
        stats!(p)
        tagsave(datadir("fig1", name), p)
    end
end

