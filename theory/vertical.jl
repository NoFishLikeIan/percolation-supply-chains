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
    using Plots, LaTeXStrings, Colors
    theme(:dao); default(size = 500 .* (√2, 1))
    plotpath = get(ENV, "PLOT_PATH", "")
end

include("definitions.jl")
include("utils/probability.jl")
include("utils/dynamicalsystem.jl")
include("utils/distcompound.jl")
include("utils/lawofmotion.jl")
include("utils/taylorapproximation.jl")

include("optimum/planner.jl")
include("optimum/agent.jl")
include("optimum/correlation.jl")

include("simulate.jl")

# Assume that m is constant and μᵢ > 0 only for i = 0.

K = 2 # Length economy
m = 100 # Size nodes
μ₀ = 0.01 # Basal risk
r = 0.01 # cost / profit ratio

model = VerticalModel(m, μ₀, r, K)

fᵣ = 0.9
sᵣ = 1.2

ρspace = 0:0.01:0.2
sspace = range(0.05, 2.5; length = 101)
vspace = range(0, m - 1; length = 101)

Eφ(ρ) = φ̃(m, sᵣ, fᵣ, ρ) / φ(m, sᵣ, 0)

plot(ρspace, Eφ)