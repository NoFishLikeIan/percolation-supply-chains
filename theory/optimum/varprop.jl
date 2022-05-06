ψ₀(x) = polygamma(0, x)
ψ₁(x) = polygamma(1, x)

G(μ, σ; sₖ, model) = G([μ, σ]; sₖ, model)
function G(x; sₖ, model::VerticalModel)
    m = model.m
    μ, σ = x
    
    pₖ = p(sₖ, μ; model = model)
    bₖ = 1 + m - μ

    ∂p = ψ₀(bₖ) - ψ₀(bₖ - sₖ)
    ∂²p = ψ₁(bₖ) - ψ₁(bₖ - sₖ)
    
    E = pₖ - (σ/2) * (1 - pₖ) * (∂p^2 + ∂²p)
    V = σ * (1 - pₖ)^2 * ∂p^2

    μ′ = m * E
    σ′ = m * (E * (1 - E) + m * V)

    return [μ′, σ′]
end

G₀(μ, m) = [
    (1 - μ) * m, (1 - μ) * μ * m 
]