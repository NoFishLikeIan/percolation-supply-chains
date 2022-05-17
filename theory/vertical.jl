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
include("utils/probability.jl")
include("utils/dynamicalsystem.jl")
include("utils/distcompound.jl")
include("utils/lawofmotion.jl")

include("optimum/planner.jl")
include("optimum/agent.jl")
include("optimum/correlation.jl")

include("simulate.jl")

# Assume that m is constant and μᵢ > 0 only for i = 0.

K = 20
m = 100
μ₀ = 0.01
profit = 100.

model = VerticalModel(m, μ₀, 1 / profit, K)
Fₛ, sₛ = compequilibrium(model)

x = [0.99, 0.4]

sspace = range(1., 2.; length = 101)
plot(sspace, s -> G₂(x; sₖ = s, model), label = nothing, c =:blue)
plot!(twinx(), sspace, s -> G₁(x; sₖ = s), label = nothing, c =:red)