function ∂ₛG₁(x; sₖ, model)
    f, ρ = x

    if ρ ≈ 0
        -(1 - f)^sₖ * log(1 - f)
    else
        (1 - model.G₁(x; sₖ)) * (
            ψ(0, sₖ + (1 - ρ) / ρ) - 
            ψ(0, sₖ + (1 - f) * (1 - ρ) / ρ)
        )
    end
end

φratio(m, s, v) = v > 0 ? φ(m, s, v) / φ(m, s, 0) : 1

function simulateφ(m, s, f, ρ; N = 15_000)
    F = BetaBinomial(m, fρtoαβ(f, ρ)...)
    F̂ = rand(F, N)
    return mean(@. φratio(m, s, F̂)) 
end

function Gfactory(m::Int64; N = 40_000)

    g(s, f, ρ) = simulateφ(m, s, f, ρ; N)

    function G₁(x; sₖ)
        f, ρ = x

        if ρ ≈ 0
            1 - (1 - f)^sₖ
        else 
            δp = (1 - ρ) / ρ
        
            1 - beta((1 - f) * δp + sₖ, f * δp) / beta((1 - f) * δp, f * δp)
        end
    end

    function G₂(x; sₖ)
        fⁿ = G₁(x; sₖ)
        Ep² = 1 - 2*(1 - fⁿ) + g(sₖ, x[1], x[2])
        return min((Ep² - fⁿ^2) / (fⁿ * (1 - fⁿ)), 1)
    end

    return G₁, G₂
end

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
