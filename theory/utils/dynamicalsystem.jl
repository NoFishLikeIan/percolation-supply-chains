"""
Check if x ∈ (0, 1) with some tolerance ε
"""
function isinopeninterval(x; ε = 1e-5)
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