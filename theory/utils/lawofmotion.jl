ψ₀(x) = polygamma(0, x)
ψ₁(x) = polygamma(1, x)

function G₁(x; s, model::VerticalModel)
    f, ρ = x
    m = model.m

    δ = ∂²p(s, m * f; model) * m * (1 - f) * f * (1 + (m - 1) * ρ)

    return p(s, m * f; model) + (δ / 2)
end

function G(x; s, model::VerticalModel)
    f, ρ = x
    m = model.m
    
    fⁿ = G₁(x; s, model)
    D = fⁿ * (1 - fⁿ)
    ρⁿ = ∂p(s, m * f; model)^2 * m * (1 - f) * f * (1 + (m - 1) * ρ) / D

    return [fⁿ, ρⁿ]
end


function sequencemoments(s::Vector{<:Real}; model::VerticalModel)
    f₀ = model.m * (1 - model.μ₀)
    ρ₀ = 0

    state = Matrix{Float64}(undef, model.K, 2)
    state[1, :] = G([f₀, ρ₀]; s = s[1], model)

    for k ∈ 2:model.K
        state[k, :] = G(state[k - 1, :]; s = s[k], model)
    end

    return state
end

JG(x; model) = ForwardDiff.jacobian(x -> G(x[1], x[2]; s = x[3], model), x)