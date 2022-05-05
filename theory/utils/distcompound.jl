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

    α, β
end
