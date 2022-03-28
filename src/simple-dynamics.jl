using Plots
using StatsBase
using Random, Distributions, Statistics
using ProgressBars
using DifferentialEquations
using LinearAlgebra


tspan = (0.0,1000)
S = 50 #number of species in the pool
α = 1 #production exponent
β = 1 #competition term
γ = 1 #other species term
d = 1 # inverse carrying capacity (diagonal element, traditional self regulation)
A = rand(Normal(1.,.1), S, S)
A[diagind(A)] .= 0 
u0 = rand(Uniform(.01,.25),S,1) 
foodweb(u,p,t) = u.^k .- u.^1.3.*A*u.^1.3.-d*u.^2.3
prob = ODEProblem(foodweb,u0,tspan)
sol = solve(prob)

plot(sol,
title = "sublinear(k=3/4),μ=1,σ=.1,K=∞,S=1,n(0)=.01",
ylabel = "n(t)",
linewidth = 2,
legend = false,
#ylims = [0,1]
)

plot(sol,
title = "logistic,μ=2,σ=1,K=1,S=1,nᵢ(0)∈[.01,.25]",
ylabel = "nᵢ(t)",
linewidth = 2,
legend = false,
ylims = [0,1]
)

savefig("Ian-multi-logistic.svg")