ψ₀(x) = polygamma(0, x)
ψ₁(x) = polygamma(1, x)

FuncDistribution = Union{Binomial, BetaBinomial}

function compoundbeta(d::Beta, m::Int64)
    α, β = params(d)
    BetaBinomial(m, α, β)
end

"""
Expected value of E[f[F]] for F ∼ BetaBinomial
"""
function E(F::FuncDistribution, f::Function)
    n = first(params(F))

    return sum( pdf(F, v) * f(v) for v ∈ 0:n )
end

"""
Variance, V[f[F]] for F ∼ BetaBinomial
"""
function V(F::FuncDistribution, f::Function)
    Ef = E(F, f)
    Eg = E(F, v -> f(v)^2)

    return Eg - Ef^2
end 

"""
Analytically matches the moment of a BetaBinomial with the μ, σ of the underlying Beta
"""
function analyticalmatchmoments(μ, σ, n)
    
    cf = -n * μ + μ^2 + σ

    α =  -cf * μ / (n*σ - n*μ + μ^2)
    β = cf * (n - μ) / (μ * (n - μ) - n*σ)

    if α < 0 || β < 0 || μ ≈ n || σ ≈ 0
        return 1., 1e-10
    else
        return α, β
    end
end

function fρtoαβ(f, ρ)
    f * (1 - ρ) / ρ, (1 - f) * (1 - ρ) / ρ
end

function αβtofρ(α, β)
    α / (α + β), 1 / (α + β + 1)
end

"""
Compute the distribution of Fₖ
"""
function inducedF(s::Real, Fₛ::FuncDistribution; model::VerticalModel)
    Epₖ = Ep(s, Fₛ; model)

    μ = model.m * Epₖ
    σ = model.m * (Epₖ * (1 - Epₖ) + model.m * Vp(s, Fₛ; model))

    α, β = analyticalmatchmoments(μ, σ, model.m)

    return BetaBinomial(model.m, α, β)

end

"""
Computes the raw moment E[F^k] of the beta-binomial distribution with (f, ρ) parameters
"""
function rawmoments(m, f, ρ)
    α, β = fρtoαβ(f, ρ)

    μ₁ = m * α / (α + β)

    μ₂ = μ₁ * (m * (1 + α) + β) / (1 + α + β)

    μ₃ = μ₁ * 
        ( m^2 * (1 + α) * (2 + α) + 3m *(1 + α) * β + β * (β - α) ) /
        (( 1 + α + β ) * ( 2 + α + β ))

    k₁ = α * m / ((α + β) * (α + β + 1) * (α + β + 2) * (α + β + 3))
    k₂ = 
        β * (α^2 - α * (4β + 1) + (β - 1) * β) +
        6 * (α + 1) * (α + 2) * β * m^2 +
        (α + 1) * β * m * (-4α + 7β - 1) +
        (α + 1) * (α + 2) * (α + 3) * m^3  

    μ₄ = k₁ * k₂

    return [μ₁, μ₂, μ₃, μ₄]
    
end