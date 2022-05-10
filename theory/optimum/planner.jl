g₁(f, σ, s; model) = G(f, σ; s, model) |> first
g₂(f, σ, s; model) = G(f, σ; s, model) |> last

function valuefunction(x::Vector{Float64}, Vₖ₊₁::Real; model::VerticalModel)
    f, σ = x
    
    J(s) = -inv(model.r) * g₁(f, σ, s; model) + s - Vₖ₊₁

    result = optimize(J, 0., 10_000.)

    sₖ = first(result.minimizer)
    Vₖ = -J(sₖ)

    sₖ, Vₖ
end

function plannerequilibrium(model::VerticalModel)

    function totalvalue(s::Vector{<:Real})
        moments = sequencemoments(s; model)

        results = Matrix{Float64}(undef, model.K, 2)
        
        for k ∈ reverse(1:model.K)
            Vₖ₊₁ = k < model.K ? results[k + 1, 2] : 0.

            results[k, :] .= valuefunction(moments[k, :], Vₖ₊₁; model)
        end

        return sum(results[:, 2])
    end

    result = optimize(s -> -totalvalue(s), ones(model.K))
    
    sₖ = result.minimizer
    moments = sequencemoments(sₖ; model)

    Fs = []

    for (μ, σ) ∈ eachrow(moments)
        α, β = analyticalmatchmoments(μ, σ, model.m)
        push!(Fs, BetaBinomial(model.m, α, β))
    end

    return Fs, sₖ
end
