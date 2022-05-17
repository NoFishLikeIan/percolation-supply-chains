ψ₀(x) = polygamma(0, x)
ψ₁(x) = polygamma(1, x)
ψ₂(x) = polygamma(2, x)

function ∂ₛG₁(x; sₖ)
    f, ρ = x

    if ρ ≈ 0
        -(1 - f)^sₖ * log(1 - f)
    else
        ForwardDiff.derivative(s -> G₁(x; sₖ = s), sₖ)
    end
end

function G₁(x; sₖ)
    f, ρ = x

    if ρ ≈ 0
        1 - (1 - f)^sₖ
    else 
        δp = (1 - ρ) / ρ
    
        1 - beta((1 - f) * δp + sₖ, f * δp) / beta((1 - f) * δp, f * δp)
    end
end

function G(x; sₖ, model::VerticalModel)
    f, ρ = x
    m = model.m

    μ = m * f
    σ = sqrt(m * f * (1 - f) * (1+ (m - 1) * ρ))

    F = BetaBinomial(m, fρtoαβ(f, ρ)...)
    γ = Distributions.skewness(F)
    κ = Distributions.kurtosis(F)
    
    fⁿ = G₁(x; sₖ)

    # https://stats.stackexchange.com/questions/452544/
    p¹ = ∂p(sₖ, μ; model)
    p² = ∂²p(sₖ, μ; model)

    a₂ = (p²^2 * μ^2 - p¹ * p² * μ + p¹^2) * σ^2
    a₃ = (-1/2) * (p¹ * p² + p²^2 * μ) * σ^3
    a₄ = (1/4) * p²^2 * σ^4

    varp = a₂ + a₃ * γ + a₄ * (κ - 1)

    ρⁿ = varp / (fⁿ * (1 - fⁿ))

    return clamp.([fⁿ, ρⁿ], 0, 1)
end

G₂(x; sₖ, model) = G(x; sₖ, model) |> last 

function sequencemoments(s::Vector{<:Real}; model::VerticalModel)
    f₀ = 1 - model.μ₀
    ρ₀ = 0.

    state = Matrix{Float64}(undef, model.K, 2)
    state[1, :] = G([f₀, ρ₀]; sₖ = s[1], model)

    for k ∈ 2:model.K
        state[k, :] = G(state[k - 1, :]; sₖ = s[k], model)
    end

    return state
end

JG(x; s, model) = ForwardDiff.jacobian(x -> G(x[1:2]; sₖ = s, model), x)[1:2, 1:2]
