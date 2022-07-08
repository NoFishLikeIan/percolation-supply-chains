using SpecialFunctions
using IterTools, Combinatorics

using Roots, Optim, NLsolve
using ForwardDiff
using LinearAlgebra
using Interpolations

using Base.Threads

using DynamicalSystems

begin
    using Random
    seed = parse(Int64, get(ENV, "SEED", "123"))
    Random.seed!(seed)

    using StatsBase, Distributions
end

begin
    using Plots, LaTeXStrings, Colors
    theme(:dao)
    default(size=500 .* (√2, 1))
    plotpath = get(ENV, "PLOT_PATH", "")
end

using DotEnv
DotEnv.config()

include("definitions.jl")

include("utils/probability.jl")
include("utils/dynamicalsystem.jl")
include("utils/distcompound.jl")
include("utils/lawofmotion.jl")
include("utils/plotting.jl")

include("optimum/planner.jl")
include("optimum/agent.jl")
include("optimum/correlation.jl")

# Assume that m is constant and μᵢ > 0 only for i = 0.

K = 20 # Length economy
m = 100 # Size nodes
μ₀ = 0.2 # Basal risk

profit = 100
costs = 3
r = costs / profit # cost / profit ratio

model = VerticalModel(m, μ₀, r, K)

X = zeros(2, model.K, 3) # (comp, soc), k, (f, ρ, s)

X[1, 1, [1, 2]] = X[2, 1, [1, 2]] = [model.μ₀, 1e-10] 

for k ∈ 2:model.K
    sagent = agentoptimum(
        X[1, k - 1, 1], X[1, k - 1, 2]; 
        m = model.m, r = model.r
    )
    splanner = planneroptimum(model, k) 


    X[1, k, [1, 2]] = G(X[1, k - 1, :]; sₖ = sagent)
    X[1, k - 1, 3] = sagent

    X[2, k, [1, 2]] = G(X[2, k - 1, :]; sₖ = splanner)
    X[2, k - 1, 3] = splanner

end

μs = X[:, :, 1]'
suppliers = X[:, :, 3]'

profits = profit * sum(1 .- μs, dims = 1) .- costs * sum(suppliers, dims = 1)


unit = range(0, 1; length = 101)

heatmap(unit, unit, (μ, ρ) -> agentoptimum(μ, ρ; m, r, integer = true), levels = 0:100)