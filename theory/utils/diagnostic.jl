function quantilerisk(q, μₖ, ρₖ; model::VerticalModel)
    if ρₖ ≈ 0 return μₖ^model.m end
    if ρₖ ≈ 1 return μₖ end

    if μₖ ≈ 1 return 1 end
    if μₖ ≈ 0 return 0 end

    α, β = fρtoαβ(1 - μₖ, ρₖ)
    F = BetaBinomial(model.m, α, β)
    x = ceil(Int64, model.m * q)

    return cdf(F, x)
end

function W̃(q, μ₀, ρ₀; model::VerticalModel)
    X = Matrix{Float64}(undef, model.K + 1, 2)
    X[1, :] = [μ₀, ρ₀]

    for k ∈ 1:model.K
        X[k + 1, :] = G̃(X[k, :], model.m, model.r)
    end

    return mean( quantilerisk(q, μ, ρ; model) for (μ, ρ) ∈ eachrow(X) )

end