struct VerticalModel
    G₁::Function
    G₂::Function
    G::Function
    m::Int64
    μ₀::Float64
    r::Float64
    K::Int64

    function VerticalModel(m, μ₀, r, K)
        G₁, G₂ = Gfactory(m)
        G(x; sₖ) = [G₁(x; sₖ), G₂(x; sₖ)]
        new(G₁, G₂, G, m, μ₀, r, K)
    end
end

