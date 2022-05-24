struct VerticalModel
    G::Function
    m::Int64
    μ₀::Float64
    r::Float64
    K::Int64

    function VerticalModel(m, μ₀, r, K)
        G = Gfactory(m)
        new(G, m, μ₀, r, K)
    end
end
