using NLsolve, Roots

using Plots, LaTeXStrings
theme(:dao)
default(size = 500 .* (√2, 1))
plotpath = "../docs/plots/"

include("definitions.jl")
include("optimum/planner.jl")
include("optimum/agent.jl")

