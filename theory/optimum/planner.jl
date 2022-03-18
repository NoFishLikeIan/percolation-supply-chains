function p(i::Int64, S::Vector{Float64})
    if i == 0 return 1. end

    risk = 1 - p(i - 1, S)
    return (1 - μ[i]) * (1 - risk)^S[i]
end

function ∂p(j, S)
    if j == 1 return 0. end

    risk = 1 - p(j - 1, S)

    return (μ[j] - 1) * log(risk) * (risk)^S[j] 
end

function D(i, j, S)
    if i == j return 1. end

    prod(
        -(1 - μ[i-k]) * S[i-k] * (1 - p(i-k-1, S))^S[i-k]
        for k ∈ 0:(i - j - 1)
    )
end

foc(Sᵣ) = foc!(Vector{Float64}(undef, length(Sᵣ)), Sᵣ)
function foc!(F, Sᵣ)
    S = [0., Sᵣ...] # Dummy for S₁
    n = length(S)
    k = length(F)

    for j ∈ 2:n
        externality = 0.
        
        for i ∈ (j + 1):n
            ratio = (m[i] * profits[i]) / (m[j] * profits[j])
            externality += ratio * D(i, j, S)
        end

        idx = j - 1
        F[idx] = ∂p(j, S) * (externality + 1) - (κ / profits[j])
    end

    return F
end