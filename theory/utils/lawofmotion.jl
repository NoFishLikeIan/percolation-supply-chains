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
    f, ρ = x

    fⁿ = G₁(x; sₖ)

    if fⁿ ≈ 1 || fⁿ ≈ 0 # By definition
        return 0.0
    end

    E₁₋ᵣ = beta((1 - ρ) / ρ, 2sₖ) / beta((1 - f) * (1 - ρ) / ρ, 2sₖ)

    return (E₁₋ᵣ - (1 - fⁿ)^2) / ((1 - fⁿ) * fⁿ)
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