using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, DataFrames
using Plots

B = [1:.1:3]
r = 1
μ = .005
K = 100
k = 3/4
Sk = 20
Sg = 3

x_eq_log(r,μ,S,K) = 1/((S-1)*μ+r/K)

x_eq_sub(r,μ,S,k) = (μ*(S-1)/r)^(1/(k-2))

# togheter

plot(B,
[[r*(1 .- x/K) for x in B]
[r*x.^(k-1) for x in B]
],
grid=false,
linewidth=3,
legends=false,
xlims = (0,100),
xticks = [0,50,100],
ylims = (0,1),
yticks = [0,.5,1],
)
plot!([[x for x in 1:.1:x_eq_log(r,μ,Sk,K)], 
[x for x in 1:.1:x_eq_log(r,μ,Sg,K)]],
[
[(Sk-1)*μ*x for x in 1:.1:x_eq_log(r,μ,Sk,K)],
[(Sg-1)*μ*x for x in 1:.1:x_eq_log(r,μ,Sg,K)]
],
grid=false,
linewidth=3,
color = :gray,
legends=false,
)
plot!(B,
[[μ*((Sk-2)*x_eq_log(r,μ,Sk,K) .+ x) for x in B]
[μ*((Sg- 2)*x_eq_log(r,μ,Sg,K) .+ x) for x in B]
[μ*((Sk-2)*x_eq_sub(r,μ,Sk,k) .+ x) for x in B]
[μ*((Sg- 2)*x_eq_sub(r,μ,Sg,k) .+ x) for x in B]
],
grid=false,
linewidth=3,
linestyle = :dash,
color = :gray,
legends=false,
)

# logistic

plot(B,
[r*(1 .- x/K) for x in B],
grid=false,
linewidth=3,
legends=false,
xlims = (0,100),
xticks = [0,50,100],
ylims = (0,1),
yticks = [0,.5,1],
)
plot!([[x for x in 1:.1:x_eq_log(r,μ,Sk,K)], 
[x for x in 1:.1:x_eq_log(r,μ,Sg,K)]],
[
[(Sk-1)*μ*x for x in 1:.1:x_eq_log(r,μ,Sk,K)],
[(Sg-1)*μ*x for x in 1:.1:x_eq_log(r,μ,Sg,K)]
],
grid=false,
linewidth=3,
color = :gray,
legends=false,
)
plot!([[x for x in x_eq_log(r,μ,Sk,K):.1:100],
[x for x in x_eq_log(r,μ,Sg,K):.1:100]],
[[μ*((Sk-2)*x_eq_log(r,μ,Sk,K) .+ x) for x in x_eq_log(r,μ,Sk,K):.1:100],
[μ*((Sg- 2)*x_eq_log(r,μ,Sg,K) .+ x) for x in x_eq_log(r,μ,Sg,K):.1:100]
],
grid=false,
linewidth=3,
linestyle = :dash,
color = :gray,
legends=false,
)
plot!([[x for x in x_eq_log(r,μ,Sk,K):.1:100],
[x for x in x_eq_log(r,μ,Sg,K):.1:100]],
[[μ*((Sk-1)*x_eq_log(r,μ,Sk,K)) for x in x_eq_log(r,μ,Sk,K):.1:100],
[μ*((Sg- 1)*x_eq_log(r,μ,Sg,K)) for x in x_eq_log(r,μ,Sg,K):.1:100]
],
grid=false,
linewidth=2,
color = :gray,
legends=false,
)

#savefig("pictorial-intra-inter-log.svg")

#vline!(1, [1],
#color = :black, 
#alpha = .5,
#linewidth=2,
#linestyle = :dash)

#plot!([x for x in 0:.1:1],
#[1 for x in 0:.1:1],
#grid=false,
#linewidth=3,
#legends=false,
#color = 2,
##xlims = (0,2),
##xticks = [0,50,100],
#ylims = (0,1.01),
##yticks = [0,.5,1],
#)

# sublinear

plot(B,
[r*x.^(k-1) for x in B],
grid=false,
linewidth=3,
legends=false,
color = 2,
xlims = (0,3),
#xticks = [0,50,100],
#ylims = (0,1),
#yticks = [0,.5,1],
)
plot!([[x for x in 1:.1:x_eq_sub(r,μ,Sk,k)], 
[x for x in 1:.1:x_eq_sub(r,μ,Sg,k)]],
[
[(Sk-1)*μ*x for x in 1:.1:x_eq_sub(r,μ,Sk,k)],
[(Sg-1)*μ*x for x in 1:.1:x_eq_sub(r,μ,Sg,k)]
],
grid=false,
linewidth=3,
color = :gray,
legends=false,
)
plot!([[x for x in x_eq_sub(r,μ,Sk,k):.1:100],
[x for x in x_eq_sub(r,μ,Sg,k):.1:100]],
[[μ*((Sk-2)*x_eq_sub(r,μ,Sk,k) .+ x) for x in x_eq_sub(r,μ,Sk,k):.1:100],
[μ*((Sg- 2)*x_eq_sub(r,μ,Sg,k) .+ x) for x in x_eq_sub(r,μ,Sg,k):.1:100]
],
grid=false,
linewidth=3,
linestyle = :dash,
color = :gray,
legends=false,
)
plot!([[x for x in x_eq_sub(r,μ,Sk,k):.1:100],
[x for x in x_eq_sub(r,μ,Sg,k):.1:100]],
[[μ*((Sk-1)*x_eq_sub(r,μ,Sk,k)) for x in x_eq_sub(r,μ,Sk,k):.1:100],
[μ*((Sg- 1)*x_eq_sub(r,μ,Sg,k)) for x in x_eq_sub(r,μ,Sg,k):.1:100]
],
grid=false,
linewidth=2,
color = :gray,
legends=false,
)
plot!([[x for x in x_eq_sub(r,μ,Sk,k):.1:100],
[x for x in x_eq_sub(r,μ,Sg,k):.1:100]],
[[(2-k)*r*x_eq_sub(r,μ,Sk,k)^(k-1)+(k-1)*r*x_eq_sub(r,μ,Sk,k)^(k-2)*x  for x in x_eq_sub(r,μ,Sk,k):.1:100],
[(2-k)*r*x_eq_sub(r,μ,Sg,k)^(k-1)+(k-1)*r*x_eq_sub(r,μ,Sg,k)^(k-2)*x for x in x_eq_sub(r,μ,Sg,k):.1:100]
],
grid=false,
linewidth=3,
linestyle = :dash,
color = 2,
legends=false,
)

savefig("pictorial-intra-inter-sub.svg")