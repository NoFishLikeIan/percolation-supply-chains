function ∂ₛG₁(x; sₖ)
    μ, ρ = x

    d = isinopeninterval(ρ) ? 
        ψ₀(sₖ + μ * (1 - ρ) / ρ) -  ψ₀(sₖ + (1 - ρ) / ρ) : 
        log(μ)

    if abs(d) < 1e-3
        return 0
    else
        return G₁(x; sₖ) * d
    end
end

"""
Agent's optimal s given the μ and ρ in previous layer.
"""
function agentoptimum(μ, ρ; m, r, integer = false)

    type = integer ? Int64 : Float64

    foc(s) = ∂ₛG₁([μ, ρ]; sₖ = s) + r

    if !isinopeninterval(μ)
        return convert(type, 1)
    end

    if foc(0) > 0 return integer ? 1 : 1e-5 end
    if foc(m) < 0 return convert(type, m) end

    s = find_zero(foc, (0, m))

    if integer
        sᵤ, sₗ = ceil(type, s), floor(type, s)

        πᵤ = G₁([μ, ρ]; sₖ = sᵤ) - r * sᵤ
        πₗ = G₁([μ, ρ]; sₖ = sₗ) - r * sₗ

        return πᵤ ≥ πₗ ? sᵤ : sₗ
    else
        return s
    end

end

"""
G with agent's sₖ = s̃
"""
function G̃(x, m::Int64, r::Float64; optkwargs...)
    sₖ = agentoptimum(x[1], x[2]; m, r, optkwargs...)
    return G(x; sₖ)
end
G̃(x, model::VerticalModel; optkwargs...) = G̃(x, model.m, model.r; optkwargs...)

function JG̃(x, m::Int64, r::Float64)
    sₖ = agentoptimum(x[1], x[2]; m, r)
    return JG(x; sₖ)
end