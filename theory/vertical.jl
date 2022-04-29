using DotEnv; DotEnv.config()

using SpecialFunctions
using IterTools, Combinatorics

using NLsolve, Roots

using Random
seed = parse(Int64, get(ENV, "SEED", "123"))
Random.seed!(seed)

using StatsBase, Distributions

using Plots, LaTeXStrings
theme(:dao); default(size = 500 .* (√2, 1))
plotpath = get(ENV, "PLOT_PATH", "")

include("definitions.jl")

include("optimum/planner.jl")
include("optimum/agent.jl")
include("optimum/correlation.jl")

include("simulate.jl")

m = 20
ratio = 0.01

F(p̃) = Binomial(m, p̃)

probspace = range(0.8, 0.99; length = 101)

function optimalsuppliers(F)
    s̃ = Vector{Float64}(undef, length(probspace))
    m̃ = convert(Float64, m)

    for (i, p̃) ∈ enumerate(probspace) 
        g(v) = pdf(F(p̃), v)
        mb(s) = ∂p(s, g; m  = m, μ = 0.)

        if mb(m̃) > ratio
            s̃[i] = m̃
        else   

            f(s) = mb(s) - ratio

            res = find_zero(f, (0., m̃))
            s̃[i] = res
        end

    end

    return s̃
end

sfig = plot(
    xlabel = L"p_k", ylabel = L"\tilde{s}",
    legendtitle = L"\kappa / \pi"
)


s̃ = optimalsuppliers(F)
label = latexstring("\$ $ratio \$")

plot!(sfig, probspace, s̃; label = label)

sfig