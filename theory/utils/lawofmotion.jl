ψ₀(x) = polygamma(0, x)
ψ₁(x) = polygamma(1, x)

G(μ, σ; s, model) = G([μ, σ]; s, model)
function G(x; s, model::VerticalModel)
    m = model.m
    μ, σ = x
    
    pₖ = p(s, μ; model = model)
    bₖ = 1 + m - μ

    ∂p = ψ₀(bₖ) - ψ₀(bₖ - s)
    ∂²p = ψ₁(bₖ) - ψ₁(bₖ - s)
    
    E = min(pₖ - (σ/2) * (1 - pₖ) * (∂p^2 + ∂²p), 1.)
    V = max(σ * (1 - pₖ)^2 * ∂p^2, 0.)

    μ′ = m * E
    σ′ = m * (E * (1 - E) + m * V)

    return [μ′, σ′]
end

G₀(μ, m) = [
    (1 - μ) * m, (1 - μ) * μ * m 
]

function sequencemoments(s::Vector{<:Real}; model::VerticalModel)
    f₀ = model.m * (1 - model.μ₀)
    σ₀ = f₀ * model.μ₀

    state = Matrix{Float64}(undef, model.K, 2)
    state[1, :] = G(f₀, σ₀; s = s[1], model)

    for k ∈ 2:model.K
        state[k, :] = G(state[k - 1, :]; s = s[k], model)
    end

    return state
end

JG(x; model) = ForwardDiff.jacobian(x -> G(x[1], x[2]; s = x[3], model), x)