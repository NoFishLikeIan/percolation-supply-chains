using DotEnv; DotEnv.config()

using SpecialFunctions
using IterTools, Combinatorics

using NLsolve, Roots

using Random
seed = parse(Int64, get(ENV, "SEED", "123"))
Random.seed!(seed)

using StatsBase, Distributions

using Plots, LaTeXStrings
theme(:dao); default(size = 500 .* (√2, 1))
plotpath = get(ENV, "PLOT_PATH", "")

include("definitions.jl")

include("optimum/planner.jl")
include("optimum/agent.jl")
include("optimum/correlation.jl")

include("simulate.jl")

function withbasalrisk(n, m, μ, profit)::VerticalModel
    VerticalModel(
        repeat([m], n), # size
        repeat([μ], n),
        repeat([profit], n),
        1.
    )
end

μ = 0.01
profit = 100
layer_size = 50
layers = 10

model = withbasalrisk(layers, layer_size, μ, profit)

G, S = solvecorrelated(model)

labels = reshape([
    latexstring("\$i = $i\$") for i ∈ 0:(layers - 1)
], (1, layers))

plot(
    0:layer_size, G';
    label = labels,
    xlabel = L"f_i", 
    ylabel = L"g_i(f_i)",
    marker = :o
)


m₀ = model.m[1]
μ₀ = model.μ[1]

F̃₀ = Normal(m₀ * (1 - μ₀), √(m₀ * (1 - μ₀) * μ₀))

plot(
    0:layer_size, G[1, :]
)

plot!(
    0:0.01:layer_size, x -> pdf(F̃₀, x)
)