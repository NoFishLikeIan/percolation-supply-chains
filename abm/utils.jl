function partitiongoods(m::Vector{Int64})

    firms = Vector{Int64}[]

    for (i, mᵢ) ∈ enumerate(m)
        baseid = i == 1 ? 1 : last(firms[i - 1]) + 1

        push!(
            firms,
            collect(baseid:(baseid - 1 + mᵢ))
        )

    end

    return firms

end

nthproduct(base, n) = Iterators.product(ntuple(i -> base, n)...)

function nthproductmatrix(base, n)
    F = Matrix{Int64}(undef, length(base)^n, n)

    for (i, tup) in enumerate(nthproduct(base, n) |> collect)
        F[i, :] .= tup
    end 

    return F

end