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
μ₀ = 0.01 # Basal risk
r = 0.1 # cost / profit ratio

model = VerticalModel(m, μ₀, r, K)

# Plots

