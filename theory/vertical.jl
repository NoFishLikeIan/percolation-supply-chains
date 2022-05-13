using DotEnv; DotEnv.config()

using SpecialFunctions
using IterTools, Combinatorics

using Roots, Optim
using ForwardDiff

begin
    using Random
    seed = parse(Int64, get(ENV, "SEED", "123"))
    Random.seed!(seed)

    using StatsBase, Distributions
end

begin
    using Plots, LaTeXStrings
    theme(:dao); default(size = 500 .* (√2, 1))
    plotpath = get(ENV, "PLOT_PATH", "")
end

include("definitions.jl")
include("utils/distcompound.jl")
include("utils/lawofmotion.jl")

include("optimum/planner.jl")
include("optimum/agent.jl")
include("optimum/correlation.jl")

include("simulate.jl")

# Assume that m is constant and μᵢ > 0 only for i = 0.

K = 3
m = 30
μ₀ = 0.01
profit = 100.

model = VerticalModel(m, μ₀, 1 / profit, K)
Fₛ, sₛ = compequilibrium(model)
Fₚ, sₚ = plannerequilibrium(model)
