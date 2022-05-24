function ∂ₛG₁(x; sₖ, model)
    f, ρ = x

    if ρ ≈ 0
        -(1 - f)^sₖ * log(1 - f)
    else
        ForwardDiff.derivative(s -> model.G(x; sₖ = s) |> first, sₖ)
    end
end

function simulateφ(m, s, f, ρ; N = 15_000)
    F = BetaBinomial(m, fρtoαβ(f, ρ)...)
    F̂ = rand(F, N)
    return mean(@. φ(m, s, F̂) / φ(m, s, 0)) 
end

function Gfactory(m::Int64; L = 15, N = 40_000)
    
    unit = range(0.01, 0.99; length = L)
    sspace = range(0.01, 5; length = L)

    paramspace = product(sspace, unit, unit)
    g(x) = simulateφ(m, x[1], x[2], x[3]; N)
    ĝ = g.(paramspace)

    φ̃ = LinearInterpolation((sspace, unit, unit), ĝ)

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
        Ep² = 1 - 2*(1 - fⁿ) + φ̃(sₖ, x[1], x[2])
        return (Ep² - fⁿ^2) / (fⁿ * (1 - fⁿ))
    end

    G(x; sₖ) = [G₁(x; sₖ), G₂(x; sₖ)]

    return G
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
