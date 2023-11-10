#= Dynamics for the consumer-resource resource framework described in the SM is defined. It is set to reproduce the results for predator-prey.
    Uncomment the relative dynamics and change parameters if necessary in order to reproduce results for consumer-resource. =#

using DrWatson, Glob
@quickactivate
include("../src/utils.jl")

using StatsBase, Random, Distributions, Statistics
using DifferentialEquations
using LinearAlgebra
using LaTeXStrings
using Plots

#= Parameters =#
nr = 50 #number of resources types / prey species 
q = 2. #consumer/resources ratio / predator/prey ratio
nc = Int64(q*nr) #number of consumer species / predator species
n = nc+nr 

# parameters for distributions
μ = 1/nr
σ = .5/√(nr)
k = μ^2/σ^2
θ = σ^2/μ
p = .5

k = 1. #exponent (1. for logistic, .75 for sublinear)

B = -ones(n); for i ∈ nc+1:n B[i] = 0 end #death vector for consumers / death vector for predators
C = ones(n); for i ∈ 1:nc C[i] = 0 end #influx vector for resources / growth vector for preys
A1, A2 = preference_matrix(n, nc; dist=Bernoulli(p)) #preference matrix (in blocks) [DEFINED AT THE BOTTOM OF THE CODE]

#= Dynamics =#
MAX_TIME = 10_000
MAX_ABUNDANCE = 1e4
blowup() = DiscreteCallback((u, t, integrator) -> maximum(u) > MAX_ABUNDANCE, terminate!)

#= community_cr(u,p,t) = u.^k .*A1*u .+ u.*A2*u .+ u .*B .+ C #dynamics for c-r model =#

community_cr(u,p,t) = u.^k .*A1*u .+ u.*A2*u .+ u .*B .+ u.^k .*C #dynamics for predator-prey model

prob = ODEProblem(
    community_cr,
    fill(.1, n), #initial condition

    (0., MAX_TIME),)

#= Solve =#
sol = solve(prob, callback = CallbackSet(TerminateSteadyState(1e-5), blowup()))

xc = [sol[i,:] for i ∈ 1:nc]
xr = [sol[i,:] for i ∈ nc+1:n]

#= Plot trajectories =#
plot(xc,
ylabel = L"x_i(t)",
xlabel = L"t",
linewidth = 2,
legend = false,
grid = false,
yscale = :log,
ylims = (1e-3,1e1),
xlims = (0,500),
palette = :Purples_9
)
plot!(xr,
linewidth = 2,
legend = false,
alpha = .33,
yscale = :log,
ylims = (1e-3,1e1),
xlims = (0,500),
palette = :Blues_9
)

#= Plot histograms =#
xc_eq = last.(xc); xc_eq[xc_eq.<=1e-4] .= 0
xr_eq = last.(xr); xr_eq[xr_eq.<=1e-4] .= 0

histogram(xc_eq, 
normalized=true,
color = "#fb8d00",
alpha = .5,
label = false,
grid = false,
#title = "Consumers distribution",
)

histogram(xr_eq, 
normalized=true,
color = :green,
alpha = .5,
label = false,
grid = false,
#title = "Resources distribution",
)

#= Preference matrix =#
function preference_matrix(n, nc; dist=(0,1))
    A1 = rand(dist, n, n) 
    A1[diagind(A1)] .= 0 
    for i ∈ 1:nc
        for j ∈ 1:nc
            A1[i,j] = 0
        end
    end
    for i ∈ nc+1:n
        for j ∈ 1:n
            A1[i,j] = 0
        end
        for j ∈ nc+1:nc
            A1[i,j] = 0
        end
    end
    A2 = zeros(n,n)
    for i ∈ nc+1:n
        for j ∈ 1:n
            A2[i,j] = -A1[j,i]
        end
    end
    return A1, A2
end