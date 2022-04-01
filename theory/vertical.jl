using NLsolve, Roots

using IterTools, Combinatorics
using StatsBase: sample

using Plots, LaTeXStrings
theme(:dao)
default(size = 500 .* (√2, 1))
plotpath = "../docs/plots/"

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
layer_size = 20
layers = 10

model = withbasalrisk(layers, layer_size, μ, profit)

G, S = solvecorrelated(m)

labels = reshape([
    latexstring("\$i = $i\$") for i ∈ 0:(layers - 1)
], (1, layers))

plot(
    0:m̄, G';
    label = labels,
    xlabel = L"f_i", 
    ylabel = L"g_i(f_i)",
    marker = :o
)
