ψ₀(x) = polygamma(0, x)
ψ₁(x) = polygamma(1, x)

G(μ, σ; s, model) = G([μ, σ]; s, model)
function G(x; s, model::VerticalModel)
    m = model.m
    μ, σ = x
    
    pₖ = p(s, μ; model = model)
    bₖ = 1 + m - μ

    ∂p = ψ₀(bₖ) - ψ₀(bₖ - s)
    ∂²p = ψ₁(bₖ) - ψ₁(bₖ - s)
    
    E = pₖ - (σ/2) * (1 - pₖ) * (∂p^2 + ∂²p)
    V = σ * (1 - pₖ)^2 * ∂p^2

    μ′ = m * E
    σ′ = m * (E * (1 - E) + m * V)

    return [μ′, σ′]
end

G₀(μ, m) = [
    (1 - μ) * m, (1 - μ) * μ * m 
]