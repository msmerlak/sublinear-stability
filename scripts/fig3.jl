using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor
using DataFrames, StatsPlots, DataFramesMeta, LaTeXStrings

@suppress_err begin 
    P = Dict{Symbol, Any}(
        :scaled => true,
        :S => 50,
        :μ => .05,
        :σ => 0.:.0001:0.003,
        :k => 1.,
        :n0 => 0.,
        :K => 100,
        :dist => "uniform",
        :N => 50,
        :seed => 1
    );


    Threads.@threads for p in expand(P)
        name = savename(p, "jld2")
        if !isfile(name)
            equilibria!(p)
            wsave(datadir("fig3", name), p)
        end
    end
end




df = collect_results(datadir("fig3"))
@df df scatter(:σ, :num_interior_equilibria)

p = Dict{Symbol, Any}(
        :scaled => false,
        :S => 500,
        :μ => .0005,
        :σ => 0.001,
        :k => .75,
        :n0 => 1.,
        :K => 1e6,
        :dist => "uniform",
        :N => 5,
        :seed => rand(UInt)
    );

evolve!(p, trajectory = true);
plot(p[:trajectory], legend = false)

boundary(p)


equilibria!(p);

p[:equilibria]
