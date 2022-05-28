function ∂ₛG₁(x; sₖ, model)
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
function agentoptimum(f, ρ; model::VerticalModel)

    bs = [0.01, model.m]

    foc(s) = ∂ₛG₁([f, ρ]; sₖ = s, model) - model.r

    roots(foc, 0..1)

    pospay = foc(bs[1]) > 0
    negmax = foc(bs[2]) < 0

    if !pospay
        bs[1]
    elseif !negmax
        bs[2]
    else
        find_zero(foc, bs)
    end
end

"""
G with agent's sₖ
"""
function G̃(x; model)
    sₖ = agentoptimum(x[1], x[2]; model)
    return G(x; sₖ)
end

function compequilibrium(model::VerticalModel)
    
    moments = Matrix{Float64}(undef, model.K, 2)
    suppliers = Vector{Float64}(undef, model.K)

    for k ∈ 1:model.K
        f, ρ = k > 1 ? moments[k - 1, :] : (1 - model.μ₀, 0.01)

        suppliers[k] = agentoptimum(f, ρ; model)
        moments[k, :] = G([f, ρ]; sₖ = suppliers[k], model)
        println("$k -> $(moments[k, :])")
    end

    return moments, suppliers
end
