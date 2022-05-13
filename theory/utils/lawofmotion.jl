ψ₀(x) = polygamma(0, x)
ψ₁(x) = polygamma(1, x)

function G₁(x; sₖ, model::VerticalModel)
    f, ρ = x
    m = model.m

    δ = ∂²p(sₖ, m * f; model) * m * (1 - f) * f * (1 + (m - 1) * ρ)

    return p(sₖ, m * f; model) + (δ / 2)
end

function G(x; sₖ, model::VerticalModel)
    f, ρ = x
    m = model.m
    
    fⁿ = G₁(x; sₖ, model)
    D = fⁿ * (1 - fⁿ)
    N = ∂p(sₖ, m * f; model)^2 * m * (1 - f) * f * (1 + (m - 1) * ρ) 

    ρⁿ = N / D

    ρⁿ = isnan(ρⁿ) ? 0 : ρⁿ

    return [fⁿ, ρⁿ]
end


function sequencemoments(s::Vector{<:Real}; model::VerticalModel)
    f₀ = 1 - model.μ₀
    ρ₀ = 0

    state = Matrix{Float64}(undef, model.K, 2)
    state[1, :] = G([f₀, ρ₀]; sₖ = s[1], model)

    for k ∈ 2:model.K
        state[k, :] = G(state[k - 1, :]; sₖ = s[k], model)
    end

    return state
end

JG(x; model) = ForwardDiff.jacobian(x -> G(x[1:2]; sₖ = x[3], model), x)