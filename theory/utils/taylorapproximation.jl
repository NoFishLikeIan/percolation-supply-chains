ψ(k, x) = polygamma(k, x)
Δψ(k, m, s, v) = ψ(k, 1 + m - v) - ψ(k, 1 + m - v - s)

function φ(m, s, v)
    if 1 + m - v - s ≤ 0
        return 0.
    end

    return (gamma(1 + m - v) / gamma(1 + m - v - s))^2
end

function a(i, m, s, v)
    if i == 1
        -2Δψ(0, m, s, v)
    elseif i == 2
        2(2Δψ(0, m, s, v)^2 - Δψ(1, m, s, v))
    elseif i == 3
        -2(2Δψ(0, m, s, v) * (2Δψ(0, m, s, v)^2 - 3Δψ(1, m, s, v)) + 
        Δψ(2, m, s, v))
        
    elseif i == 4

        big = 1 + m - v
        small = 1 + m - v - s

        ΔH = Δψ(0, m, s, v)

        c = 8ΔH^4 - 24ψ(1, big) * ΔH^2 + 8ψ(2, big) * ΔH +
        6ψ(1, small) * (4ΔH^2 + ψ(1, small) - 2ψ(1, big)) -
        8ΔH * ψ(2, small) + 
        ψ(3, small) + 6ψ(1, big)^2 - ψ(3, big)
        
        return 2c
    end
end

"""
Taylor approximation of E[φratio] around 0
"""
function φratio(m, s, f, ρ)

    μₛ = rawmoments(m, f, ρ)

    1 + sum(a(i, m, s, 0) * μₛ[i] / factorial(i) for i ∈ 1:4)
    
end