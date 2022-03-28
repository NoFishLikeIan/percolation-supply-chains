using NLsolve, Roots

using IterTools
using StatsBase: sample

using Plots, LaTeXStrings
theme(:dao)
default(size = 500 .* (√2, 1))
plotpath = "../docs/plots/"

include("definitions.jl")
include("optimum/planner.jl")
include("optimum/agent.jl")
include("simulate.jl")

n = 3
	
m = VerticalModel(
    repeat([100], n), # m
    [0. for i ∈ 1:n], # μ
    100. .* ones(n), # π
    1 # κ
)

foc! = focfactory(m)
foc(S) = foc!(Vector{Float64}(undef, n - 1), S)


res = nlsolve(foc!, ones(n-1))
@assert res.f_converged

Sₛ = [1., res.zero...]
Z, _ = integersolution(Sₛ, m)

F = resilience(ones(Int64, n); m, T = 10_000)