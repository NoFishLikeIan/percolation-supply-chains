function plannerequilibrium(model::VerticalModel)
    r, μ₀ = model.r, model.μ₀
    m, K = model.m, model.K

    if r > (1 - μ₀^m)
        return zeros(K)
    else
        sₖ = ones(K)
        sₖ[1] = m
        return sₖ
    end
end

function plannersymmetricsolution(model::VerticalModel; L=20, verbose=false)
    m, r = model.m, model.r
    fspace = range(0.01, 0.99; length=L)
    ρspace = copy(fspace)

    sspace = range(0, m; length=L)

    sₖ = Array{Float64}(undef, (L, L, model.K))
    Vₖ = copy(sₖ)

    indices = product(1:L, 1:L) |> collect

    for k ∈ reverse(1:K)
        verbose && print("Planner: solving layer $k / $K...\r")

        Threads.@threads for (i, j) ∈ indices
            f, ρ = fspace[i], ρspace[j]

            if k < K

                Jₖᵢⱼ = Vector{Float64}(undef, length(sspace))
                for (l, s) ∈ enumerate(sspace)
                    f′, ρ′ = G([f, ρ]; sₖ=s)
                    itp = LinearInterpolation(
                        (fspace, ρspace), Vₖ[:, :, k+1];
                        extrapolation_bc=Line()
                    )

                    Jₖᵢⱼ[l] = model.m * (f′ - s * model.r) + itp(f′, ρ′)
                end

                idxmax = argmax(Jₖᵢⱼ)
                sₖ[i, j, k] = sspace[idxmax]
                Vₖ[i, j, k] = Jₖᵢⱼ[idxmax]

            else
                s = agentoptimum(f, ρ; m, r)
                f′ = G₁([f, ρ]; sₖ=s)

                sₖ[i, j, end] = s
                Vₖ[i, j, end] = model.m * (inv(model.r) * f′ - s)

            end
        end
    end

    function planneroptimum(f, ρ, k::Int64)
        LinearInterpolation(
            (fspace, ρspace), sₖ[:, :, k];
            extrapolation_bc=Line()
        )(f, ρ)
    end

    return planneroptimum
end