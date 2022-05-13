ψ₀(x) = polygamma(0, x)
ψ₁(x) = polygamma(1, x)

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
    
    fⁿ = G₁(x; sₖ)
    D = fⁿ * (1 - fⁿ)
    N = ∂p(sₖ, m * f; model)^2 * m * (1 - f) * f * (1 + (m - 1) * ρ) 

    ρⁿ = N / D

    ρⁿ = isnan(ρⁿ) ? 0 : ρⁿ

    return [fⁿ, ρⁿ]
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