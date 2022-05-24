function ∂ₛG₁(x; sₖ)
    f, ρ = x

    if ρ ≈ 0
        -(1 - f)^sₖ * log(1 - f)
    else
        ForwardDiff.derivative(s -> G₁(x; sₖ = s), sₖ)
    end
end

function simulateφ(m, s, f, ρ; N = 15_000)
    F = BetaBinomial(m, fρtoαβ(f, ρ)...)
    F̂ = rand(F, N)
    return mean(@. φ(m, s, F̂) / φ(m, s, 0)) 
end

function Gfactory(
    model::VerticalModel; 
    order = 10, 
    lb = [0.01, 0.01, 0.01],
    ub = [5., 0.99, 0.99], N = 40_000)

    paramspace = chebpoints((order, order, order), lb, ub)

    g(x) = simulateφ(model.m, x[1], x[2], x[3]; N)
    φ̃ = chebinterp(g.(paramspace), lb, ub)

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
        Ep² = 1 - 2*(1 - fⁿ) + φ̃([sₖ, x[1], x[2]])
        return (Ep² - fⁿ^2) / (fⁿ * (1 - fⁿ))
    end

    G(x; sₖ) = [G₁(x; sₖ), G₂(x; sₖ)]

    return G
end

function sequencemoments(s::Vector{<:Real}; model::VerticalModel)
    f₀ = 1 - model.μ₀
    ρ₀ = 0.
    G = Gfactory(model)

    state = Matrix{Float64}(undef, model.K, 2)
    state[1, :] = G([f₀, ρ₀]; sₖ = s[1], model)

    for k ∈ 2:model.K
        state[k, :] = G(state[k - 1, :]; sₖ = s[k], model)
    end

    return state
end

JG(x; s, model) = ForwardDiff.jacobian(x -> G(x[1:2]; sₖ = s, model), x)[1:2, 1:2]
