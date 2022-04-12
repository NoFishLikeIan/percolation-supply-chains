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

function allsolutions!(Z::Matrix{Int64}, model::VerticalModel)

    n = length(model.m)
	
	# Social planner
	foc! = focfactory(model)
	res = nlsolve(foc!, ones(n - 1))

	if res.f_converged
		Sₛ = [1, res.zero...] 
		Zₛ, _ = integersolution(Sₛ, model)
        Z[:, 1] = Zₛ
    else
        throw("Convergence error for μ = $(model.μ[1]) and r = $(model.profits[1] |> inv)")
	end

	# Competitive w/out correlation

	_, Zₐ = sequence(n; m = model, integer = true)
	Zₐ = collect(Int64, Zₐ)
    Z[:, 2] = Zₐ

	# Competitive w/ correlation
	_, Zᵨ = solvecorrelated(model)
    Z[:, 3] = Zᵨ	
end


function paramphase(μs::Vector{Float64}, rs::Vector{Float64}; verbose = true)

    Z = ones(Int64, layers, 3)
    M, N = length(μs), length(rs)

    cached = Dict{Vector{Int64}, Float64}()
    R = Array{Float64}(undef, M, N, 3)

    model = withbasalrisk(layers, layer_size, μ, profit)

    for i ∈ 1:M
        model.μ[1] = μs[i]
        
        for j ∈ 1:N    
            verbose && print("$(M*(i - 1) + j) / $(M * N) \r")

            model.profits .= inv(rs[j])
            
            allsolutions!(Z, model)

            for t ∈ 1:3
                Zₜ = Z[:, t]
                cachedres = get(cached, Zₜ, nothing)

                if isnothing(cachedres)
                    res = resilience(Zₜ; m = model, T = 5_000) |> mean
                    cached[Zₜ] = res

                    R[i, j, t] = res
                else
                    R[i, j, t] = cachedres
                end
            end
        end
    end

    return R

end
