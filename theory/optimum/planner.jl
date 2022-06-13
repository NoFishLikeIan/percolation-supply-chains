function plannerequilibrium(model::VerticalModel; L = 101, verbose = false)
    m, r = model.m, model.r
    fspace = range(0.01, 0.99; length = L)
    ρspace = copy(fspace)

    sspace = range(0.01, m - 1; length = 2L)
    
    sₖ = Array{Float64}(undef, (L, L, model.K))
    Vₖ = Array{Float64}(undef, (L, L, model.K))

    for i ∈ 1:L, j ∈ 1:L
        f, ρ = fspace[i], ρspace[j]
        s = agentoptimum(f, ρ; m, r)
        f′ = G₁([f, ρ]; sₖ = s)

        sₖ[i, j, end] = s
        Vₖ[i, j, end] = model.m * ( inv(model.r) * f′ - s )
    end
    
    for k ∈ reverse(1:(K-1))
        verbose && print("Layer $k...\r")
        
        for i ∈ 1:L, j ∈ 1:L

            f, ρ = fspace[i], ρspace[j]
            Jₖᵢⱼ = Vector{Float64}(undef, 2L)  
            for (l, s) ∈ enumerate(sspace) 
                f′, ρ′ = G([f, ρ]; sₖ = s)
                itp = LinearInterpolation(
                    (fspace, ρspace), Vₖ[:, :, k + 1]; 
                    extrapolation_bc = Line()
                )

                Jₖᵢⱼ[l] = model.m * ( inv(model.r) * f′ - s ) + itp(f′, ρ′)
            end
            
            idxmax = argmax(Jₖᵢⱼ)
            sₖ[i, j, k] = sspace[idxmax]
            Vₖ[i, j, k] = Jₖᵢⱼ[idxmax]
        end
    end

    return sₖ, Vₖ
end
