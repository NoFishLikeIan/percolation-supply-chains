ψ₀(x) = polygamma(0, x)
ψ₁(x) = polygamma(1, x)
ψ₂(x) = polygamma(2, x)


function G₁(x; sₖ)
    f, ρ = x

    if ρ ≈ 0
        1 - (1 - f)^sₖ
    else 
        δρ = (1 - ρ) / ρ
    
        1 - beta((1 - f) * δρ + sₖ, f * δρ) / beta((1 - f) * δρ, f * δρ)
    end
end

function G(x; sₖ, model::VerticalModel)
    f, ρ = x
    m = model.m

    VF = m * (1 - f) * f * (1 + (m - 1) * ρ)
    SF = VF * (1 - 2f) * (1 + (2m - 1) * ρ) / (1 + ρ)

	NKF = ρ^2 * (6*(m-1)*m + 2 - 3f*(7*(m-1)*m + 2)*(1-f)) + 3*(1-f)*f*(((m-8)*m + 4)*ρ + m - 2) + 6m*ρ - 3ρ + 1
	DKF = (ρ + 1) * (2ρ + 1)
	KF = VF * NKF / DKF
    
    fⁿ = G₁(x; sₖ)

    order2 = ∂p(sₖ, m * f; model)^2 * VF
    order3 = ∂p(sₖ, m * f; model) * ∂²p(sₖ, m * f; model) * SF
    order4 = (∂²p(sₖ, m * f; model)^2 / 4 + (2 / 9) * ∂p(sₖ, m * f; model) * ∂³p(sₖ, m * f; model)) * KF

    ρⁿ = (order2 + order3 + order4) / (fⁿ * (1 - fⁿ))

    return clamp.([fⁿ, ρⁿ], 0, 1)
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

JG(x; model) = ForwardDiff.jacobian(x -> G(x[1:2]; sₖ = x[3], model), x)
