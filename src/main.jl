# Packages
using Agents
using Graphs
using StatsBase

# Plotting
using Plots, GraphPlot

default(size = 500 .* (√2, 1))
theme(:dao)

# Random
using Random

rng = Random.MersenneTwister(1212)

# Subfiles
include("utils.jl")
include("definitions.jl")
include("agentstep.jl")

include("plots.jl")

m = [3, 3, 3]
μ = 0.05 * ones(sum(m))

function verticaleconomy(m::Vector{Int64}, μ::Vector{Float64})

    n = sum(m)
    goods = partitiongoods(m)

    model = ABM(
        Firm, 
        properties = Dict(:goods => goods),
        scheduler = Schedulers.by_id
    )

    for i ∈ 1:n 
        firm = Firm(i, [], true, μ[i])
        s = potentialsuppliers(firm, goods)

        firm.x = isempty(s) ? Int64[] : [sample(s)]
        add_agent!(firm, model) 
    end

    return model

end

model = verticaleconomy(m, μ)
T = 30

data, _ = run!(model, agent_step!, dummystep, T; adata = [:isfunctional, :x])


function plott(t)

    tindex = data.step .== t
    functional = data[tindex, :isfunctional]
    x = data[tindex, :x]
    plotgraph(partitiongoods(m), functional, x)

end
