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
function Sr(i, pₛ; m::VerticalModel, integer = false, brackets = [1., 100.])
    ratio = m.κ / ((1 - m.μ[i]) * m.profits[i])
    risk = 1 - pₛ

    f = S -> log(risk) * (risk)^S + ratio

    lb, ub = brackets
    isbracketing = f(lb) * f(ub) < 0

    S̄ = isbracketing ? find_zero(f, brackets) : find_zero(f, 2.)

    if !integer return S̄ end

    prof(s) = m.profits[i] * p(i, pₛ, s; m) - m.κ * s 

    if prof(ceil(S̄)) > prof(floor(S̄))
        return ceil(S̄)
    else
        return floor(S̄)
    end

end

function p(i::Int64, pₛ::Float64, s; m::VerticalModel)
    if i == 1 return (1 - m.μ[i]) end

    supres = 1 - pₛ
    return (1 - m.μ[i]) * (1 - supres^s)
end

function p(i::Int64, pₛ::Float64; m::VerticalModel, integer = false)
    if i == 1 return (1 - m.μ[i]) end
    s = Sr(i, pₛ; m, integer)
    return p(i, pₛ, s; m)
end


function sequence(L::Int64; m::VerticalModel, integer = false)
    S = ones(L)
    ps = Vector{Float64}(undef, L); ps[1] = p(1, 1.; m, integer)

    for l ∈ 2:L 
        ps[l] = p(l, ps[l - 1]; m, integer)
        S[l] = Sr(l, ps[l - 1]; m, integer)
    end

    return ps, S
end