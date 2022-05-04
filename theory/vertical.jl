using DotEnv; DotEnv.config()

using SpecialFunctions
using IterTools, Combinatorics
using DynamicalSystems

using NLsolve, Roots

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

include("optimum/planner.jl")
include("optimum/agent.jl")
include("optimum/correlation.jl")
include("optimum/varprop.jl")

include("simulate.jl")

# Assume that m is constant and μᵢ > 0 only for i = 0.

K = m = 20
μs = zeros(K); μs[1] = 0.01
profit = 100.

model = VerticalModel(m, 0.01, 0.01)


function makeds(model, s)

    function varprop!(dx, x, params, k)
        s = first(params)

        sₖ = s # Maybe -> k > 1 ? 1. : s

        dx .= G(x; s = sₖ, model = model, k = k + 1)
    end

    return DiscreteDynamicalSystem(varprop!, x₀, [s])
end

function getfinal(μ₀, s; K = 20)

    model = VerticalModel(
        repeat([m], K), # m
        [μ₀, zeros(K - 1)...],
        repeat([profit], K), 1.
    )

    ds = makeds(model, s)
	
	tr = trajectory(ds, K)

	notnan = findfirst(row -> any(isnan.(row)), eachrow(tr))

	idx = K + 1 #isnothing(notnan) ? K + 1 : notnan - 1

	μ, σ = tr[idx, :]
	
	return μ, max(0, σ)
end

S = 401
μspace = range(0.01, 0.7; length = S)
sspace = range(0.1, 2.; length = S)
