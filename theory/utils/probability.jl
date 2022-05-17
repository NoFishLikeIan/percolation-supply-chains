
"""
Using Sterling approximation of factorial
"""
function p̃(s::Real, v::Real; model::VerticalModel)
    m = model.m
    Δₛ(x) = x * log((x - s) / x) - s * log(x - s)

    1 - exp(Δₛ(m) - Δₛ(m - v))
end

function p(s::Real, v::Real; model::VerticalModel)
    m = model.m
    
    nodiv = s == 0
    fulldiv = s > 1 + m - v
    largem = m > 100

    if nodiv return 0. end
    if fulldiv return 1. end
    if largem return p̃(s, v; model) end

    constant = gamma(1 + m - s) / gamma(1 + m)
    random = gamma(1 + m - s - v)  / gamma(1 + m - v)

    return 1. - constant / random    
end

function ∂p(s::Real, v::Real; model::VerticalModel)
    m = model.m

    Δψ₀ = ψ₀(1 + m - v) - ψ₀(1 + m - v - s)
    return (1 - p(s, v; model)) * Δψ₀
end

function ∂²p(s::Real, v::Real; model::VerticalModel)
    m = model.m

    Δψ₀ = ψ₀(1 + m - v) - ψ₀(1 + m - v - s)
    Δψ₁ = ψ₁(1 + m - v) - ψ₁(1 + m - v - s)
    
    return -(1 - p(s, v; model)) * (Δψ₀^2 + Δψ₁)
end

function ∂³p(s::Real, v::Real; model::VerticalModel)
    m = model.m

    Δψ₀ = ψ₀(1 + m - v) - ψ₀(1 + m - v - s)
    Δψ₁ = ψ₁(1 + m - v) - ψ₁(1 + m - v - s)
    Δψ₂ = ψ₂(1 + m - v) - ψ₂(1 + m - v - s)

    pₖ = p(s, v; model)
    pₖ′ = ∂p(s, v; model)
    

    return pₖ′ * (Δψ₀^2 + Δψ₁) + (1 - pₖ) * (2Δψ₀ * Δψ₁ + Δψ₂)
end