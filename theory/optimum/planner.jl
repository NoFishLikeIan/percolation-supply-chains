function p(i::Int64, S::Vector{Float64}; m::VerticalModel)
    if i == 0 return 1. end

    risk = 1 - p(i - 1, S; m)
    return (1 - m.μ[i]) * (1 - risk)^S[i]
end

function ∂p(i, j, S; m::VerticalModel)
    if i < j || j == 1 return 0. end

    ownrisk = 1 - m.μ[i]
    uprisk = 1 - p(i - 1, S; m)

    if i == j
        return -ownrisk * uprisk^S[i] * log(uprisk)
    elseif i > j
        return ownrisk * S[i] * uprisk^S[i - 1] * ∂p(i - 1, j, S; m)
    end
end

function E(i, j, S; m::VerticalModel)
    if i ≤ j return 0. end
    if i == j + 1 return 1. end

    uprisk = 1 - p(i-1, S; m)

    return S[i] * (1 - m.μ[i]) * uprisk^(S[i] - 1) * E(i-1, j, S; m)
end

function focfactory(m::VerticalModel)
    function foc!(F, Sᵣ)
        S = [0., Sᵣ...] # Dummy for S₀

        n = length(S)

        for j ∈ 2:n

            ext = sum(
                E(i, j, S; m) * (m.m[i] * m.profits[i]) / (m.m[j] * m.profits[j])
                for i ∈ 1:n
            )

            F[j - 1] = ∂p(j, j, S; m) * (ext + 1) - (m.κ / m.profits[j])
        end

        return F
    end
end