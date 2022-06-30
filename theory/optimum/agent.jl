function ∂ₛG₁(x; sₖ)
    f, ρ = x

    if ρ ≈ 0
        -(1 - f)^sₖ * log(1 - f)
    else
        (1 - G₁(x; sₖ)) * (
            ψ(0, sₖ + (1 - ρ) / ρ) -
            ψ(0, sₖ + (1 - f) * (1 - ρ) / ρ)
        )
    end
end

"""
Agent's optimal s given the f and ρ in previous layer.
"""
function agentoptimum(f, ρ; m, r)

    brackets = [0.01, m]

    foc(s) = ∂ₛG₁([f, ρ]; sₖ=s) - r

    pospay = foc(brackets[1]) > 0
    negmax = foc(brackets[2]) < 0

    if !pospay
        brackets[1]
    elseif !negmax
        brackets[2]
    else
        find_zero(foc, brackets)
    end
end

"""
G with agent's sₖ
"""
function G̃(x, m::Int64, r::Float64)
    sₖ = agentoptimum(x[1], x[2]; m, r)
    return G(x; sₖ)
end
G̃(x, model::VerticalModel) = G̃(x, model.m, model.r)

function JG̃(x, m::Int64, r::Float64)
    sₖ = agentoptimum(x[1], x[2]; m, r)
    return JG(x; sₖ)
end