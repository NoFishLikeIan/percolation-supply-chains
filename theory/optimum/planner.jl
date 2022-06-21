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
