

function agent_step!(firm::Firm, model::AgentBasedModel)
    s = potentialsuppliers(firm, model.goods)

    isbasal = isempty(s)

    ispotfunctional = isbasal || any( model.agents[i].isfunctional for i ∈ firm.x )

    firm.isfunctional = ispotfunctional && rand(rng) > firm.μ
end