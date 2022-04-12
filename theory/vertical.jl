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

layer_size = 50
layers = 15

parameters = 100

μs = range(1e-3, 0.05; length = parameters) |> collect
rs = range(0.001, 0.1; length = parameters) |> collect

R = paramphase(μs, rs)
ΔR = R[:, :, 1] .- R[:, :, 2]

cmax = maximum(abs, extrema(ΔR))

contourf(
    μs, rs, ΔR';
    xlims = extrema(μs), ylims = extrema(rs),
    clims = (-cmax, cmax),
    xlabel = L"\mu", ylabel = L"\kappa / \pi",
    title = "Social planner",
)