Goods = Vector{Vector{Int64}}

mutable struct Firm <: AbstractAgent
    id::Int64

    θ::Vector{Float64} # Perceived probability over neighbors
    μ::Float64 # Idiosyncratic risk

    isfunctional::Bool
    x::Vector{Bool} # Boolean array of suppliers
    λ::Float64 # Learning rate

    p::Float64
    k::Float64

    realized::Float64 # Realized profits
end

function good(f::Firm, G::Goods)
    findfirst(g -> f.id ∈ g, G)
end

function potentialsuppliers(f::Firm, goods::Goods)
    g = good(f, goods)

    return get(goods, g - 1, Int64[])
end
