function p(i::Int64, s::Int64, f::Int64; m::VerticalModel)
    idx = i + 1
    mᵢ = m.m[idx]

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
    
    function nextprob(f)
            
        total = 0.

        for fₚ ∈ 0:mₚ
            prob = upstreamprob[fₚ + 1]
            total += gₚ[fₚ + 1] * prob^f * (1 - prob)^(mᵢ - f)
        end

        return binomial(mᵢ, f) * total
    end

    return nextprob.(fs)

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