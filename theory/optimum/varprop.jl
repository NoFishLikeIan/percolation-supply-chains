ψ₀(x) = polygamma(0, x)
ψ₁(x) = polygamma(1, x)

function p̃(k , s, v; model::VerticalModel)

    idx = k # Just for clarity, since k - 1 = idx but indexing in julia starts at 1.
    m = model.m[idx]

    if s > 1 + m - v return 1. end

    constant = gamma(1 + m - s) / gamma(1 + m)
    random = gamma(1 + m - v - s) / gamma(1 + m - v)

    return 1 - constant / random
end

G(μ, σ; s, model, k) = G([μ, σ]; s, model, k)
function G(
    x; 
    s, model::VerticalModel, k::Int64
)
    idx = k
    m = model.m[idx]

    μ, σ = x
    
    pₖ = p̃(k, s, μ; model = model)
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