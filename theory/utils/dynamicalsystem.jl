"""
f nullcline for G̃
"""
function NG̃f(ρ, model::VerticalModel)
    find_zeros(f -> first(G̃([f, ρ], model.m, model.r)) - f, 0, 1)
end

"""
ρ nullcline for G̃
"""
function NG̃ρ(f, model::VerticalModel)
    find_zeros(ρ -> last(G̃([f, ρ], model.m, model.r)) - ρ, 0, 1)
end

function fixedpoints(model, x₀)
    function f!(F, x)
        F .= G̃([f, ρ], model.m, model.r) - x
    end
    
    sol = mcpsolve(
        f!, [0., 0.], [1., 1.], x₀,
        autodiff = :forward
    )

    return sol.zero

end 