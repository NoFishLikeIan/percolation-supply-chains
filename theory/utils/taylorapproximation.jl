ψ(k, x) = polygamma(k, x)
Δψ(k, m, s, v) = ψ(k, 1 + m - v) - ψ(k, 1 + m - v - s)

φ(m, s, v) = (gamma(1 + m - v) / gamma(1 + m - v - s))^2
φ¹(m, s, v) = -2φ(m, s, v) * Δψ(0, m, s, v)
φ²(m, s, v) = 2φ(m, s, v) * (2Δψ(0, m, s, v)^2 + Δψ(1, m, s, v))
φ³(m, s, v) = -2φ(m, s, v) * (4Δψ(0, m, s, v)^3 + 6Δψ(0, m, s, v) * Δψ(1, m, s, v) + Δψ(2, m, s, v))

function φ⁴(m, s, v)
    u = 1 + m - v
    l = 1 + m - v - s
    
    coeff = 8Δψ(0, m, s, v)^4 + 
        24ψ(1, u)*Δψ(0, m, s, v)^2 +
        8ψ(2, u)*Δψ(0, m, s, v) + 
        6ψ(1, l)*(
            -4Δψ(0, m, s, v)^2 + ψ(1, l) - 2ψ(1, u)
        ) -
        8ψ(2, l)*Δψ(0, m, s, v) - 
        ψ(3, l) + 6ψ(1, u)^2 + ψ(3, l)

    return 2φ(m, s, v) * coeff
end

"""
Taylor approximation of the φ function
"""
function φₜ(m, s, v₀, Δv₁, Δv₂, Δv₃, Δv₄)	
    φ(m, s, v₀) + 
    φ¹(m, s, v₀) * Δv₁ + 
    φ²(m, s, v₀) * Δv₂^2 / 2 + 
    φ³(m, s, v₀) * Δv₃^3 / 6 + 
    φ⁴(m, s, v₀) * Δv₄^4 / 24
end
function φₜ(m, s, v, v₀)
    Δv = (v - v₀)
    φₜ(m, s, v₀, Δv, Δv, Δv, Δv)
end

"""
Taylor approximation of φ with around E[F]
"""
function φ̃(m, s, f, ρ)

    μ, σ², γ, κ = moments(m, f, ρ)
    σ = √σ²

    return φₜ(
        m, s, 
        μ, 0, 
        σ² * σ, γ * σ^3, κ * σ^4
    )
end