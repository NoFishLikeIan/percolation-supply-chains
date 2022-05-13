"""
Agent's optimal s given the choice of s in the previous node and the ratio of cost and expected profits.
"""
function sₒ(sₛ::Real, Fₛ::FuncDistribution; model::VerticalModel, integer = false, br = [0.01, 100.])

    function foc(s)
        function ∂pₛ(v) # FIXME: Is this the correct expectation?
            pₛ = p(sₛ, v; model)
            if pₛ < 1.
                (1 - pₛ)^s * log(1 - pₛ)    
            else
                ∂pₛ(v - 1) # Not sure...
            end       
        end

        return E(Fₛ, ∂pₛ) + model.r
    end

    lb, ub = br

    s̄ = foc(ub) * foc(lb) < 0 ? find_zero(foc, br) : find_zero(foc, 0.05)

    if !integer 
        return s̄
    end

    # Assumes unit costs
    prof(s) = inv(model.r) * E(Fₛ, v -> p(s, v; model)) - s

    if prof(ceil(s̄)) > prof(floor(s̄))
        return ceil(Int64, s̄)
    else
        return floor(Int64, s̄)
    end

end

function Ep(s::Real, Fₛ::FuncDistribution; model::VerticalModel)
    f(v) = p(s, v; model)
    return E(Fₛ, f)
end
function Vp(s::Real, Fₛ::FuncDistribution; model::VerticalModel)
    f(v) = p(s, v; model)
    return V(Fₛ, f)
end


function p(s::Real, v::Real; model::VerticalModel)
    if s == 0
        return 0.
    end

    m = model.m

    if s > 1 + m - v
        return 1.
    else
        constant = gamma(1 + m - s) / gamma(1 + m)
        random = gamma(1 + m - s - v)  / gamma(1 + m - v)

        return 1. - constant / random
    end    
end

function ∂p(s::Real, v::Real; model::VerticalModel)
    m = model.m
    Δψ₀ = ψ₀(1 + m - v) - ψ₀(1 + m - v - s)
    return (1 - p(s, v; model)) * Δψ₀
end

function ∂²p(s::Real, v::Real; model::VerticalModel)
    m = model.m

    Δψ₀ = ψ₀(1 + m - v) - ψ₀(1 + m - v - s)
    Δψ₁ = ψ₁(1 + m - v) - ψ₁(1 + m - v - s)
    return -(1 - p(s, v; model)) * (Δψ₀^2 + Δψ₁)
end


function compequilibrium(model::VerticalModel; integer = false)
    F₀ = Binomial(model.m, 1 - model.μ₀)

    # No correlation in the 0th, hence 1st admits an analytical solution
    firstsup = (log(model.r) - log(-log(model.μ₀))) / log(model.μ₀)

    s₁ = integer ? round(Int64, firstsup) : firstsup
    F₁ = inducedF(s₁, F₀; model)

    Fs = Vector{FuncDistribution}(undef, model.K); Fs[1] = F₁
    suppliers = Vector{integer ? Int64 : Float64}(undef, model.K)
    suppliers[1] = s₁

    for k ∈ 2:model.K
        Fₛ = Fs[k - 1]
        sₛ = suppliers[k - 1]

        suppliers[k] = sₒ(sₛ, Fₛ; model, integer)
        Fs[k] = inducedF(suppliers[k], Fₛ; model)
    end

    moments = [ αβtofρ(params(F)[2:3]...) for F ∈ Fs ]

    return moments, suppliers
end
