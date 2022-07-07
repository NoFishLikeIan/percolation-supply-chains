function planneroptimum(m::Int64, r::Float64, μ₀::Float64, k::Int64)
    n = (log(r) - (log ∘ log)(1/μ₀)) / log(μ₀)

    if n < 0 return 0 end
    if k > 1 return 1 end

    return min(n, m)
end

"""
G with planner's sₖ = s̃ₚ
"""
function Gₚ(x, m::Int64, r::Float64, μ₀::Float64, k::Int64)
    sₖ = planneroptimum(m, r, μ₀, k)
    return G(x; sₖ)
end
Gₚ(x, model::VerticalModel, k::Int64) = Gₚ(x, model.m, model.r, model.μ₀, k)

function JGₚ(x, m::Int64, r::Float64, μ₀::Float64, k::Int64)
    sₖ = planneroptimum(m, r, μ₀, k)
    return JG(x; sₖ)
end