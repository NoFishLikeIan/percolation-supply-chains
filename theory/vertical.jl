using NLsolve, Roots

using Plots, LaTeXStrings
theme(:dao)
default(size = 500 .* (√2, 1))

include("optimum/agent.jl")
include("optimum/planner.jl")

m = [200 for i ∈ 1:9]
n = length(m)
μ = 0.01 .* ones(n)
profits = 100. .* ones(n)
κ = 1.

# Social planner vs agent, S
res = nlsolve(foc!, ones(n - 1))
@assert res.f_converged
Ssocial = res.zero

pagent, Sagent = sequence(n)
compfig = plot(xlabel = "Layer", ylabel = L"S_i")

plot!(compfig, 2:n, Ssocial; 
    marker = :o, label = "Social planner")

plot!(compfig, 2:n, Sagent[2:n]; 
    marker = :o, label = "Firm problem")

compfig