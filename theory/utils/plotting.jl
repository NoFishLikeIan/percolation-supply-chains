function plotvectorfield(xs, ys, g::Function, plotargs...; plotkwargs...)
    fig = plot()
    plotvectorfield!(fig, xs, ys, g, plotargs...; plotkwargs...)
    return fig
end

function plotvectorfield!(figure, xs, ys, g::Function, plotargs...; rescale = 1, plotkwargs...)

    N, M = length(xs), length(ys)
    Δx = last(extrema(xs)) - first(extrema(xs))

    xm, ym = (repeat(xs, outer=length(ys)), repeat(ys, inner=length(xs)))

    field = g.(xm, ym)

    scale = rescale * Δx / min(N, M)
    u = @. scale * first(field)
    v = @. scale * last(field)

    steadystates = @. (u ≈ 0) * (v ≈ 0)

    u[steadystates] .= NaN
    v[steadystates] .= NaN

    Δx = maximum(xs) - minimum(xs)
    Δy = maximum(ys) - minimum(ys)

    xlims = (minimum(xs) - Δx / 20, maximum(xs) + Δx / 20)
    ylims = (minimum(ys) - Δy / 20, maximum(ys) + Δy / 20)

    quiver!(
        figure, xm, ym, plotargs...;
        quiver = (u, v), line_z=repeat(norm.(field), inner=4),
        aspect_ratio = 1, xlims = xlims, ylims = ylims,
        c = :batlow,
        plotkwargs...
    )

end

function agentvectorfieldplot(m, r; L = 20, rescale = 2)
    Φ(μ, ρ) = G̃([μ, ρ], m, r) .- [μ, ρ]
    μspace = range(0, 1; length = L)
    ρspace = range(0, 1 - 1e-3; length = L)


    fig = plotvectorfield(
        μspace, ρspace, Φ; 
        xlabel = L"\mu", ylabel = L"\rho", 
        title = latexstring("Vector field \$ \\tilde{G}(x) - x \$, with \$\\kappa / \\pi = $(r)\$"), rescale = rescale
    )

    plot!(fig,
        range(extrema(μspace)...; length = 501),
        μ -> steadystatecurve(μ, r), c = :darkred,
        label = nothing
    )

    scatter!(
        fig, find_zeros(x -> x*log(x) + r, (0, 1)), zeros(2);
        label = nothing, c = [:darkred, :darkblue]
    )

    scatter!(
        fig, [0, 1], zeros(2);
        label = nothing, c = [:darkblue, :darkred]
    )


    return fig
end