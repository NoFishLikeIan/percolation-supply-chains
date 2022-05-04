ψ₀(x) = polygamma(0, x)
ψ₁(x) = polygamma(1, x)

FuncDistribution = Union{Binomial, BetaBinomial}

function compoundbeta(d::Beta, m::Int64)
    α, β = params(d)
    BetaBinomial(m, α, β)
end

"""
Expected value of E[f[F]] for F ∼ BetaBinomial
"""
function E(F::FuncDistribution, f::Function)
    n = first(params(F))

    return sum( pdf(F, v) * f(v) for v ∈ 0:n )
end

"""
Variance, V[f[F]] for F ∼ BetaBinomial
"""
function V(F::FuncDistribution, f::Function)
    Ef = E(F, f)
    Eg = E(F, v -> f(v)^2)

    return Eg - Ef^2
end 

"""
Compute the distribution of Fₖ
"""
function inducedF(s::Real, Fₛ::FuncDistribution; model::VerticalModel)
    m = model.m
    Epₖ = Ep(s, Fₛ; model)

    μ = m * Epₖ
    σ = m * (Epₖ * (1 - Epₖ) + m * Vp(s, Fₛ; model))

    function momentsmatch!(F, p)
        α, β = p

        F[1] = μ - m * α / (α + β) 
        F[2] = σ - m * α * β * (α + β + m) / ((α + β)^2 * (α + β + 1))
    end

    res = nlsolve(momentsmatch!, [0.5, 0.5])
    α, β = res.zero

    return BetaBinomial(m, α, β)

end


"""
Agent's optimal s given the choice of s in the previous node and the ratio of cost and expected profits.
"""
function sₒ(sₛ::Real, Fₛ::FuncDistribution; model::VerticalModel, integer = false, br = [0.01, 100.])

    function foc(s)
        function g(v) # FIXME: Is this the correct expectation?
            pₛ = p(sₛ, v; model)
            (1 - pₛ)^s * log(1 - pₛ)            
        end

        return s * E(Fₛ, g) - model.r
    end

    lb, ub = br
    isbracketing = foc(lb) * foc(ub) < 0

    s̄ = isbracketing ? find_zero(foc, br) : find_zero(foc, 2.)

    if !integer 
        return s̄
    end

    # Assumes unit costs
    p̄(v) = p(s̄, v; model)
    prof(s) = inv(model.r) * E(p, Fₛ) - s

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


function p(s::Real, v::Integer; model::VerticalModel)
    m = model.m

    if s > 1 + m - v
        return 1.
    else
        constant = gamma(1 + m - s) / gamma(1 + m)
        random = gamma(1 + m - s - v)  / gamma(1 + m - v)

        return 1. - constant / random
    end    
end

function ∂p(s::Real, v::Integer; model::VerticalModel)
    m = model.m
    if s > 1 + m - v
        return 0.
    else
        (1 - p(s, v; model)) * (ψ₀(1 + m - v) - ψ₀(1 + m - v - s))
    end 
end


function equilibrium(K::Int64; model::VerticalModel, integer = false)
    F₀ = Binomial(model.m, 1 - model.μ)

    s₁ = (log(model.r) - log(-log(model.μ₀))) / log(model.μ₀)
    F₁ = inducedF(s₁, F₀; model)

    Fs = Vector{FuncDistribution}(undef, K); Fs[1] = F₁
    suppliers = ones(K); suppliers[1] = s₁

    for k ∈ 2:(K + 1)
        Fₛ = Fs[k - 1]
        sₛ = suppliers[k - 1]

        suppliers[k] = sₒ(sₛ, Fₛ; model, integer)
        Fs[k] = inducedF(sₛ, Fₛ; model)

    end

    return probs, suppliers
end