function model_step!(model)

    print("Step $(model.t) \r")
    model.t += 1

    for firms ∈ model.goods, i ∈ firms        
        picksuppliers!(model[i], model)
        updatebeliefs!(model[i], model)
        drawfunctional!(model[i], model)
        computepayoffs!(model[i], model)
    end
end

function getF(firm::Firm, model)
    s = potentialsuppliers(firm, model.goods)
    return [ model[i].isfunctional ? 1 : 0 for i in s ]
end

function computepayoffs!(firm::Firm, model)
    F = getF(firm, model)
    firm.realized = firm.p * max.(F'firm.x) - firm.k * sum(firm.x)
end

function updatebeliefs!(firm::Firm, model)
    if isempty(firm.x) return end

    F = getF(firm, model)
    g = good(firm, model.goods)

    s = findfirst(
        row -> all(row .== F), 
        eachrow(model.Fₛ[g - 1]) |> collect
    )

    firm.θ[s] = (firm.θ[s] + 1) / 2

    firm.θ[1:s-1] = firm.θ[1:s-1] ./ 2
    firm.θ[s+1:end] = firm.θ[s+1:end] ./ 2

end

function picksuppliers!(firm::Firm, model)
    if isempty(firm.x) return end

    g = good(firm, model.goods)
    X = model.Fₛ[g - 1]

    Π = Vector{Float64}(undef, size(X, 1))

    for (i, x) ∈ eachrow(X) |> enumerate
        Π[i] =  firm.p * max.(X*x, 1)'firm.θ - firm.k * sum(x)
    end

    firm.x = X[argmax(Π), :]
    
end

function drawfunctional!(firm::Firm, model)

    F = getF(firm, model)

    ispotfunctional = !isempty(F) ? F'firm.x > 0 : true

    firm.isfunctional = ispotfunctional && rand(model.rng) > firm.μ
end