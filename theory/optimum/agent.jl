"""
Agents optimal's S over the integers given the probability of the previous' node and the ratio of cost and expected profits.
"""
function Sz(pₛ, profit, κ, μ)
    ratio = κ / ((1 - μ) * profit)
    return (log(ratio) - log(pₛ)) / log(1 - pₛ) 
end

"""
Agents optimal's S over the reals given the probability of the previous' node and the ratio of cost and expected profits.
"""
function Sr(pₛ, profit, κ, μ)
    ratio = κ / ((1 - μ) * profit)
    risk = 1 - pₛ

    f = S -> log(risk) * (risk)^S + ratio

    return find_zero(f, [1, 10])
end

function p(i::Int64, pₛ::Float64)
    if i == 1 return (1 - μ[i]) end

    Sᵢ = Sr(pₛ, profits[i], κ, μ[i])
    risk = (1 - pₛ)

    return (1 - μ[i]) * (1 - risk)^Sᵢ
end

function sequence(L::Int64)
    S = ones(L)
    ps = Vector{Float64}(undef, L); ps[1] = p(1, 1.)

    for l ∈ 2:L 
        ps[l] = p(l, ps[l - 1])
        S[l] = Sr(ps[l - 1], profits[l], κ, μ[l])
    end

    return ps, S
end