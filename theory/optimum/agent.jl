function ∂ₛG₁(x; sₖ)
    μ, ρ = x

    d = isinopeninterval(ρ) ? 
        ψ₀(sₖ + μ * (1 - ρ) / ρ) -  ψ₀(sₖ + (1 - ρ) / ρ) : 
        log(μ)

    return G₁(x; sₖ) * d
end

"""
Agent's optimal s given the μ and ρ in previous layer.
"""
function agentoptimum(μ, ρ; m, r)

    foc(s) = ∂ₛG₁([μ, ρ]; sₖ = s) + r

    try
        if foc(0) > 0 || !isinopeninterval(μ)
            1e-5
        elseif foc(m) < 0
            m
        else
            find_zero(foc, (0, m))
        end
    catch e
        println("μ, ρ = $μ, $ρ")

        throw(e)
    end
end

"""
G with agent's sₖ = s̃
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