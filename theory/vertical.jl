using NLsolve, Roots

using Plots, LaTeXStrings
theme(:dao)
default(size = 500 .* (√2, 1))
plotpath = "../docs/plots/"

include("definitions.jl")
include("optimum/planner.jl")
include("optimum/agent.jl")

n = 40
	
m = VerticalModel(
    repeat([10], n), # m
    [0.01 for i ∈ 1:n], # μ
    100. .* ones(n), # π
    1 # κ
)

foc! = focfactory(m)
foc(S) = foc!(Vector{Float64}(undef, n - 1), S)


res = nlsolve(foc!, ones(n-1))
@assert res.f_converged