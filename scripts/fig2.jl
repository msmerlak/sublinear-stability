using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor

@suppress_err begin 
    P = Dict{Symbol, Any}(
        :scaled => false,
        :S => 50,
        :μ => .001:.001:0.01,
        :σ => .001:.001:0.01,
        :k => .75,
        :n0 => 1.,
        :K => 1000.,
        :dist => "uniform",
        :N => 10,
        :seed => 1
    );

    Threads.@threads for p in expand(P)
        name = savename(p, "jld2")
        stats!(p)
        tagsave(datadir("fig2", name), p)
    end
end

