using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using DataFrames, StatsPlots, DataFramesMeta, LaTeXStrings


df = collect_results(datadir("fig2"))

@df df scatter(
    :μ * 50,
    :σ,
    marker_z = :richness,
    xlabel = "μ"
)
