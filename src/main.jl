# Packages
using Agents
using Graphs
using StatsBase
using Combinatorics
using IterTools

# Plotting
using Plots, GraphPlot

default(size = 500 .* (√2, 1))
theme(:dao)

# Random
using Random

# Subfiles
include("utils.jl")
include("definitions.jl")
include("step.jl")

include("plots.jl")

m = [3, 3, 3] # Market sizes
μ = 0.05 * ones(sum(m)) # Idyosincratic risk    
θ = [ones(2^mᵢ) ./ 2^mᵢ for mᵢ ∈ m[2:end]] # Uniform priors
Fₛ = [
    nthproductmatrix([0, 1], mᵢ) for mᵢ ∈ m[2:end]
]

p, k = 10., 1.

function verticaleconomy(
    m::Vector{Int64}, μ::Vector{Float64};
    seed = 1212
)

    rng = Random.MersenneTwister(seed)
    goods = partitiongoods(m)

    properties = Dict(
        :goods => goods,
        :Fₛ => Fₛ, # State space
        :rng => rng
    )

    model = ABM(Firm; properties = properties)

    for (g, good) ∈ enumerate(goods)
        isbasal = g == 1

        θ₀ = isbasal ? Float64[] : θ[g - 1]
        x₀ = isbasal ? Int64[] : repeat([true], m[g - 1])

        for i ∈ good
            firm = Firm(i, θ₀, μ[i], true, x₀, p, k)
            add_agent!(firm, model) 
        end

    end

    return model

end

model = verticaleconomy(m, μ)

T = 1_000

data, _ = run!(model, dummystep, model_step!, T; adata = [:isfunctional, :x, :θ])


if false


    function plott(t)

        tindex = data.step .== t
        functional = data[tindex, :isfunctional]
        x = data[tindex, :x]
        plotgraph(partitiongoods(m), functional, x)

    end
end

