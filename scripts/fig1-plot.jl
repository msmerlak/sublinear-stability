using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using DataFrames, StatsPlots, DataFramesMeta, LaTeXStrings

df = collect_results(datadir("fig1"))


X = @subset(df, :n0 .== 1.)
sort!(X, :S)
plot(dpi = 500, title = "with threshold")
for K in [20, 50, 100, 500]
    @df @subset(X, (:k .== 1) .& (:K .== K)) plot!(
        :S, 
        :richness, 
        yerror = :diversity_se,
        dpi = 500,
        color = COLOR_LOG,
        label = "K = $K",
        markershape = :auto,
        xlabel = "species richness " * L"S",
        ylabel = "fraction of surviving species "* L"\Phi"
        )
end

@df @subset(X, (:k .== .75) .& (:K .== 500)) plot!(
    :S, 
    :richness,
    color = COLOR_SUB, 
    yerror = :diversity_se,
    markershape = :auto,
    label = "sublinear"
    )

current()
savefig(plotsdir("fig1-threshold.png"))


X = @subset(df, :n0 .== 0.)
sort!(X, :S)
plot(dpi = 500, title = "without threshold")
for K in [20, 50, 100, 500]
    @df @subset(X, (:k .== 1) .& (:K .== K)) plot!(
        :S, 
        :richness, 
        yerror = :diversity_se,
        dpi = 500,
        color = COLOR_LOG,
        label = "K = $K",
        markershape = :auto,
        xlabel = "species richness " * L"S",
        ylabel = "fraction of surviving species "* L"\Phi"
        )
end

@df @subset(X, (:k .== .75) .& (:K .== 500)) plot!(
    :S, 
    :richness,
    color = COLOR_SUB, 
    yerror = :diversity_se,
    markershape = :auto,
    label = "sublinear"
    )

current()
savefig(plotsdir("fig1-threshold.png"))
