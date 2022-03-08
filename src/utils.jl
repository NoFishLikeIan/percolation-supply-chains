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