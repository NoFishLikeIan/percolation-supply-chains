Goods = Vector{Vector{Int64}}

mutable struct Firm <: AbstractAgent
    id::Int64
    isfunctional::Bool # i ∈ F

    x::Vector{Int64} # Suppliers
    b::Matrix{Float64} # Perceived probability over neighbors
    μ::Float64 # Idiosyncratic risk
end

function good(f::Firm, G::Goods)
    findfirst(g -> f.id ∈ g, G)
end

function potentialsuppliers(f::Firm, goods::Goods)
    g = good(f, goods)

    return get(goods, g - 1, Int64[])
end
