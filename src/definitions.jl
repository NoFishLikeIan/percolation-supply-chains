Goods = Vector{Vector{Int64}}

mutable struct Firm <: AbstractAgent
    id::Int64
    x::Vector{Int64}
    isfunctional::Bool
    μ::Float64
end

function good(f::Firm, G::Goods)
    findfirst(g -> f.id ∈ g, G)
end

function potentialsuppliers(f::Firm, goods::Goods)
    g = good(f, goods)

    return get(goods, g - 1, Int64[])
end
