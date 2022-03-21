"""
Agents optimal's S over the integers given the probability of the previous' node and the ratio of cost and expected profits.
"""
function Sz(i, pₛ; m::VerticalModel)
    ratio = m.κ / ((1 - m.μ[i]) * m.profits[i])
    return (log(ratio) - log(pₛ)) / log(1 - pₛ) 
end

"""
Agents optimal's S over the reals given the probability of the previous' node and the ratio of cost and expected profits.
"""
function Sr(i, pₛ; m::VerticalModel)
    ratio = m.κ / ((1 - m.μ[i]) * m.profits[i])
    risk = 1 - pₛ

    f = S -> log(risk) * (risk)^S + ratio

    return find_zero(f, [1., 10.])
end

function p(i::Int64, pₛ::Float64; m::VerticalModel)
    if i == 1 return (1 - m.μ[i]) end

    Sᵢ = Sr(i, pₛ; m)
    risk = (1 - pₛ)

    return (1 - m.μ[i]) * (1 - risk)^Sᵢ
end

function sequence(L::Int64; m::VerticalModel)
    S = ones(L)
    ps = Vector{Float64}(undef, L); ps[1] = p(1, 1.; m)

    for l ∈ 2:L 
        ps[l] = p(l, ps[l - 1]; m)
        S[l] = Sr(l, ps[l - 1]; m)
    end

    return ps, S
end