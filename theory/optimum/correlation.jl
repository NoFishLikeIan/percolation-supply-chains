function p(i::Int64, s::Int64, f::Int64; m::VerticalModel)
    idx = i + 1
    mᵢ = m.m[idx - 1]

    return (1 - m.μ[idx]) * (1 - binomial(mᵢ-s, f) / binomial(mᵢ, f))
end

function correlateds(gₚ::Vector{Float64}, i::Int64; m::VerticalModel)
    idx = i + 1
    mₚ = m.m[idx-1]

    for s ∈ 0:mₚ
        prob = sum(gₚ[f + 1] * (p(i, s + 1, f; m) - p(i, s, f; m)) for f ∈ 0:mₚ)

        if prob * m.profits[idx] < m.κ
            return s
        end
    end

    return mⱼ
end

"""
Computes the next distribution g, assuming s
"""
function gᵢ(i::Int64, gₚ::Vector{Float64}, s::Int64; m::VerticalModel)
    idx = i + 1

    mᵢ = m.m[idx]
    mₚ = m.m[idx-1]
    fs = 0:mᵢ

    upstreamprob = p.(i, s, 0:mₚ; m)
    
    function expectedprob(f)
            
        E = 0.

        for fₚ ∈ 0:mₚ
            prob = upstreamprob[fₚ + 1]
            E += gₚ[fₚ + 1] * prob^f * (1 - prob)^(mᵢ - f)
        end

        return binomial(mᵢ, f) * E
    end

    return expectedprob.(fs)

end

function gᵢ(i::Int64, gₚ::Vector{Float64}; m::VerticalModel)
    s = correlateds(gₚ, i; m)

    return gᵢ(i, gₚ, s; m)
end

function solvecorrelated(m::VerticalModel)
    m̄ = maximum(m.m)
    layers = length(m.m)

    G = Matrix{Float64}(undef, (layers, m̄ + 1))
    μ₀ = m.μ[1]

    G[1, :] .= [ 
        binomial(m.m[1], f) * (1 - μ₀)^f * μ₀^(m.m[1] - f) 
        for f ∈ 0:m.m[1]         
    ]

    S = [1]

    for l ∈ 2:layers
        i = l - 1
        s = correlateds(G[l - 1, :], i; m)
        G[l, :] = gᵢ(i, G[l - 1, :], s; m)

        push!(S, s)
    end

    return G, S
end

function ∂p(i::Int64, s::Float64, f::Int64; m::VerticalModel)
    idx = i + 1
    mₚ = m.m[idx - 1]
    μᵢ = m.μ[idx]

    num = gamma(1 + mₚ - s) * gamma(1 + mₚ - f) 
    den = gamma(1 + mₚ - s - f) * gamma(1 + mₚ)

    H = log(mₚ - s) - log(mₚ - f - s)

    return (1 - μᵢ) * H * num / den

end