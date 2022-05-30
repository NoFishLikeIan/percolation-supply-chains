"""
f nullcline for G̃
"""
function NG̃f(ρ; model)
    find_zeros(f -> G̃([f, ρ]; model)[1] - f, 0, 1)
end

"""
ρ nullcline for G̃
"""
function NG̃ρ(f; model)
    find_zeros(ρ -> G̃([f, ρ]; model)[2] - ρ, 0, 1)
end

function fixedpoints(model, x₀)
    function f!(F, x)
        F .= G̃(x; model) - x
    end
    
    sol = mcpsolve(
        f!, [0., 0.], [1., 1.], x₀,
        autodiff = :forward
    )

    return sol.zero

end 