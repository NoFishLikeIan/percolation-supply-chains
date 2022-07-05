"""
Check if x ∈ (0, 1)
"""
function isinopeninterval(x; ε = 0)
    ε < x < 1 - ε
end

function steadystatecurve(μ, r)
    function focss(ρ)
        der = ∂ₛG₁([μ, ρ]; sₖ = 1)
        return der + r
    end

    try
        find_zero(focss, (1e-3, 1 - 1e-3))

    catch
        NaN
    end
end