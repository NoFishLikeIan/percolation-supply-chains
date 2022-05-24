"""
Agent's optimal s given the f and ρ in previous layer.
"""
function agentoptimum(f, ρ; model::VerticalModel)

    foc(s) = ∂ₛG₁([f, ρ]; sₖ = s, model) - model.r

    sₖ = foc(model.m) < 0 ? find_zero(foc, [0.01, model.m]) : model.m

    return sₖ
end

"""
G with agent's sₖ
"""
function G̃(x; model)
    sₖ = agentoptimum(x[1], x[2]; model)
    return model.G(x; sₖ)
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
