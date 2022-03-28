function resilience(Z::Vector{Int64}; m::VerticalModel, T = 10_000)

    N = maximum(m.m)
    L = length(m.m)

    F = Array{Int64}(undef, T, L, N)

    for t ∈ 1:T

        F[t, 1, :] = collect(Int64, rand(first(m.m)) .> first(m.μ))

        for l ∈ 2:L, i ∈ 1:m.m[l]

            suppliers = sample(1:m.m[l - 1], Z[l])
            PF = sum(F[t, l - 1, suppliers]) > 0
            
            F[t, l, i] = (PF && rand() > m.μ[l]) ? 1 : 0

        end
    end


    return F
end