# Packages
using Agents
using Graphs
using StatsBase
using Combinatorics
using IterTools
using Random

# Plotting
using Plots, GraphPlot

default(size = 500 .* (√2, 1))
theme(:dao)


# Utilities
include("utils.jl")

# Agents
include("definitions.jl")
include("step.jl")

# Plotting
include("plots.jl")

m = [10, 5, 5, 5] # Market sizes
μ = 0.01 * ones(sum(m)) # Idyosincratic risk    

θ = [ ones(2^mᵢ) ./ 2^mᵢ for mᵢ ∈ m[1:end-1] ] # Uniform priors
Fₛ = [ nthproductmatrix([0, 1], mᵢ) for mᵢ ∈ m[1:end - 1] ]
payoff = ones(sum(m))

function verticaleconomy(
    m::Vector{Int64}, μ::Vector{Float64}, payoff::Vector{Float64},
    θ::Vector{Vector{Float64}}, Fₛ::Vector{Matrix{Int64}};
    k = 1, seed = 1212
)

    rng = Random.MersenneTwister(seed)
    goods = partitiongoods(m)

    properties = Dict(
        :goods => goods, :Fₛ => Fₛ, # State space
        :t => 1, :rng => rng
    )

    model = ABM(Firm; properties = properties)

    for (g, good) ∈ enumerate(goods)
        isbasal = g == 1

        θ₀ = isbasal ? Float64[] : θ[g - 1]
        x₀ = isbasal ? Int64[] : repeat([true], m[g - 1])

        for i ∈ good
            firm = Firm(i, θ₀, μ[i], true, x₀, payoff[i], k, 0.)
            add_agent!(firm, model) 
        end

    end

    return model

end

model = verticaleconomy(m, μ, payoff, θ, Fₛ)

T = 250

data, _ = run!(
    model, dummystep, model_step!, T; 
    adata = [:isfunctional, :x, :θ]
)