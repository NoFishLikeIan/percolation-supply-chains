function G₁(x; sₖ)
    μ, ρ = x

    if ρ ≈ 0
        return μ^sₖ
    end

    ρᵣ = (1 - ρ) / ρ
    return beta(μ * ρᵣ + sₖ, (1 - μ) * ρᵣ) / beta(μ * ρᵣ, (1 - μ) * ρᵣ)
end

function G₂(x; sₖ)
    μ, ρ = x
    μ′ = G₁(x; sₖ)

    if μ′ ≈ 1 || μ′ ≈ 0 || ρ ≈ 0 # By definition
        return 0.
    end

    ρᵣ = (1 - ρ) / ρ
    E₁₋ᵣ = beta(ρᵣ, 2sₖ) / beta(μ * ρᵣ, 2sₖ)
    
    return (E₁₋ᵣ - (μ′)^2) / (μ′ * (1 - μ′))
end

G(x; sₖ) = [G₁(x; sₖ), G₂(x; sₖ)]

function JG(x; sₖ)

    return ForwardDiff.jacobian(v -> G(v; sₖ), x)

    # TODO: Finish analytical implementation
    if false
        x′ = G(x; sₖ)

        f, ρ = x
        f′, ρ′ = x′
        α, β = fρtoαβ(f, ρ)

        ρ̃ = (1 - ρ) / ρ

        G₁₁ = (1 - f′) * ρ̃ * (ψ(0, β + sₖ) - ψ(0, β))
        G₁₂ = (1 - f′) *
              (
            (1 - f) * (ψ(0, β + sₖ) - ψ(0, β)) -
            (ψ(0, ρ̃ + sₖ) - ψ(0, ρ̃))
        ) / ρ^2

        N = (ρ′ + (1 - f′) / f′) * f′ * (1 - f′)
        ∂N = -N * ρ̃ * (ψ(0, β + 2sₖ) - ψ(0, β))
    end

end

function sequencemoments(s::Vector{<:Real}; model::VerticalModel)
    f₀ = 1 - model.μ₀
    ρ₀ = 0.0

    state = Matrix{Float64}(undef, model.K, 2)
    state[1, :] = G([f₀, ρ₀]; sₖ=s[1], model)

    for k ∈ 2:model.K
        state[k, :] = G(state[k-1, :]; sₖ=s[k], model)
    end

    return state
end