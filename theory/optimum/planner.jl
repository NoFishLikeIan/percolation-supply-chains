function p(i::Int64, S::Vector; m::VerticalModel)
    if i == 0 return 1. end

    supres = 1 - p(i - 1, S; m)
    ownrisk = 1 - m.μ[i]

    return ownrisk * (1 - supres^S[i])
end

function ∂p(i, j, S; m::VerticalModel)
    if i < j || j == 1 return 0. end

    ownrisk = 1 - m.μ[i]
    supres = 1 - p(i - 1, S; m)

    if i == j
        return -ownrisk * supres^S[i] * log(supres)
    elseif i > j
        return ownrisk * S[i] * supres^S[i - 1] * ∂p(i - 1, j, S; m)
    end
end

function E(i, j, S; m::VerticalModel)
    if i ≤ j return 0. end
    if i == j + 1 return 1. end

    uprisk = 1 - p(i-1, S; m)

    return S[i] * (1 - m.μ[i]) * uprisk^(S[i] - 1) * E(i-1, j, S; m)
end

function focfactory(m::VerticalModel)
    function foc!(F, Sᵣ)
        S = [1., Sᵣ...] # Dummy for S₀

        n = length(S)

        for j ∈ 2:n

            ext = sum(
                E(i, j, S; m) * (m.m[i] * m.profits[i]) / (m.m[j] * m.profits[j])
                for i ∈ 1:n
            )

            F[j - 1] = ∂p(j, j, S; m) * (ext + 1) - (m.κ / m.profits[j])
        end

        return F
    end
end

"""
Given a solution over R, computes the solution over Z
"""
function integersolution(S::Vector{Float64}, m::VerticalModel)
    n = length(S)

    if n > 20
        throw("n=$n too big")
    end

    function profit(Z)
        probs = (i -> p(i, Z; m)).(1:n)
        return @. m.profits * probs - m.κ * Z
    end

    space = product(((floor(Int64, s), ceil(Int64, s)) for s ∈ S)...)

    maxvec = Vector{Int64}(undef, n)
    maxprof = -Inf

    for Ztup ∈ space

        Z = collect(Ztup)
        πₛ = profit(Z) |> sum
        if πₛ > maxprof
            maxprof = πₛ
            maxvec .= Z
        end
    end

    return maxvec, maxprof

end
