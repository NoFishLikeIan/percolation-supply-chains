### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8b32f51d-e1f7-489d-95af-5e25909d709d
using PlutoUI

# ╔═╡ 963896ed-ae03-4b68-bec2-6b9917093454
# ╠═╡ show_logs = false
begin
	using Plots, LaTeXStrings, Colors
    theme(:dao)

	FONT_SIZE = 14
		
	default(
		size = 650 .* (√2, 1),
		legend = :topright, dpi = 300, fmt = :eps,
		legendfontsize = FONT_SIZE, legendtitlefontsize = FONT_SIZE,
		xtickfontsize = FONT_SIZE, ytickfontsize = FONT_SIZE,
		xguidefontsize = FONT_SIZE, yguidefontsize = FONT_SIZE,
		linewidth = 2.5
	)
end

# ╔═╡ edf0c6f6-db46-11ec-19d0-f3a901cbdf5e
# ╠═╡ show_logs = false
begin
	using SpecialFunctions
	using IterTools, Combinatorics
	using LinearAlgebra
	
	using Roots, Optim, NLsolve
	using ForwardDiff
	using Interpolations

	using Random
    using Distributions
	
    Random.seed!(123)
    
end

# ╔═╡ 9f11ee9e-e89f-41c4-8ef4-d91c7d4e8db3
# ╠═╡ show_logs = false
using DynamicalSystems

# ╔═╡ c87845ef-2a91-41c8-933a-b6e139739927
begin
	include("../theory/definitions.jl")
	include("../theory/utils/probability.jl")
	include("../theory/utils/dynamicalsystem.jl")
	include("../theory/utils/distcompound.jl")
	include("../theory/utils/lawofmotion.jl")
	include("../theory/utils/plotting.jl")
	include("../theory/utils/diagnostic.jl")
	
	include("../theory/optimum/planner.jl")
	include("../theory/optimum/agent.jl")
	include("../theory/optimum/correlation.jl")
	
end

# ╔═╡ d3add18c-4845-4683-a82a-cbbe16f32b6f
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ 1a39488c-43c8-47e3-a760-913962cb3703
SAVE = false

# ╔═╡ 82d52cc7-c6d4-45e5-8eeb-e81fda6317f8
function makecurvefrompoints(points)
	T = length(points)
	C = maximum(length.(points))

	curve = NaN .* zeros(T, C)

	for (t, p) ∈ enumerate(points)
		for (i, pᵢ) ∈ enumerate(p)
			curve[t, i] = pᵢ
		end
	end

	return curve
end

# ╔═╡ 551ba67b-5317-4dd0-9b9d-126ca8c13eb0
m = 100; K = 50

# ╔═╡ 077b426a-ac05-4b6d-bcf1-0e379a82ba76
Norbit, Ttr = 2000, 2000

# ╔═╡ 7bb7dd41-0768-476f-b815-536607ecbdf0
unit = range(0.01, 0.99; length = 50)

# ╔═╡ 7253f132-6b20-43aa-8df8-fd8ef15cbcf5
md"
## Interpretation of Beta distribution
"

# ╔═╡ a0ff3a4e-c189-4459-b610-1a7109ffbcbc
begin

	function betacdf(μ, ρ)
		α, β = fρtoαβ(1 - μ, ρ)
		d = Beta(α, β)

		return x -> cdf(d, x)
	end

	function betabinpdf(μ, ρ)
		α, β = fρtoαβ(1 - μ, ρ)
		d = BetaBinomial(1_000, α, β)

		return x -> pdf(d, x)
	end
end

# ╔═╡ ab5e7adf-8acc-4fff-92f5-c5d66d23f7f7
begin
	betadomain = 0:0.001:1
	μfig = 0.5
	ρs = [0.001, 0.1, 0.5] 

	intensity_colors = [:darkred, :darkorange, :darkblue]
	
	betafig = plot(
		xlabel = L"p", ylabel = latexstring("Cumulative density, \$\\mathbb{P} (P \\leq p)\$"),
		margins = 5Plots.mm, legend = :topleft, legendtitle = L"\rho"
	)

	for (i, ρ) ∈ ρs |> enumerate
		plot!(betafig,
			betadomain, betacdf(μfig, ρ); 
			c = intensity_colors[i], label = latexstring("\$$(ρ)\$")
		)
	end

	vline!(
		betafig, [1 - μfig],
		linestyle = :dash, c = :black, label = nothing
	)

	annotate!(
		betafig, 
		1 - μfig, 0.1, text(L"\ \mathbb{E}[P] = 1 - \mu", :black, 18, :left)
	)
	

	

	betafig

end

# ╔═╡ 7c2dafe1-13ce-4c57-8c83-ef0d4b235302
begin
	binomialdomain = 0:1:1_000
	
	binomialfig = plot(
		;
		xlabel = L"n", ylabel = latexstring("Probability density, \$\\mathbb{P} (F = n)\$"),
		margins = 5Plots.mm, legend = :topleft, legendtitle = L"\rho", label = L"0"
	)

	for (i, ρ) ∈ ρs |> enumerate
		plot!(binomialfig,
			binomialdomain, betabinpdf(μfig, ρ); 
			c = intensity_colors[i], label = latexstring("\$$(ρ)\$")
		)
	end
	
	binomialfig

end

# ╔═╡ e8397cfc-407a-451d-b571-d95dd66c8c33
if SAVE 
	savefig(betafig, joinpath("../docs/plots", "beta-cdf"))
	savefig(binomialfig, joinpath("../docs/plots", "betabin-pdf"))
end

# ╔═╡ 15f7b2fb-f98c-464d-9055-59bb796adb7b
md"
## Size of dumpening
"

# ╔═╡ 1e41fda9-6e2e-4402-b685-72d49e6caf9a
begin
	dumpfig = plot(
		xlabel = latexstring("\$n\$-th additional supplier"), 
		ylabel = "Risk dumpening factor",
		margins = 4Plots.mm, legend = :right,
		legendtitle = L"\rho"
	)

	for (i, ρ) ∈ enumerate(ρs)
		
		plot!(dumpfig,
			0:1:20,
			marker = :o,
			n -> inv(1 + n * ρ / (1 - ρ)),
			label = latexstring("\$$(ρ)\$"),
			c = intensity_colors[i]
		)

	end


	dumpfig
end

# ╔═╡ 5c1f36ac-2d44-4011-af8c-7634473ffc54
if SAVE savefig(dumpfig, joinpath("../docs/plots", "risk-dumpening")) end

# ╔═╡ b67eb0ba-6db7-40ba-a0b2-21151c4eb748
md"
## Function $G$
"

# ╔═╡ 49e1c995-831b-472e-b68c-b56a21d1f308
md"
``s`` $(@bind sₖ Slider(
	0.5:0.01:4, show_value = true, default = 1.3
))
"

# ╔═╡ e518339b-9781-4f37-895f-19343315c79c
begin
	L = 13
	μspace = range(1e-3, 0.99; length = L)
	ρspace = range(1e-3, 0.99; length = L)
	space = product(μspace, ρspace)
end

# ╔═╡ 5228005b-19f6-4cab-bd31-31aedb4fcf6c
let

	dμ = (x -> first(G([x[1], x[2]]; sₖ)) - x[1]).(space)
	dρ = (x -> last(G([x[1], x[2]]; sₖ)) - x[2]).(space)

	μlim = maximum(abs.(dμ))
		
	g1fig = contourf(
		μspace, ρspace, dμ';
		c = :coolwarm, xlabel = L"\mu", ylabel = L"\rho", title = L"G_{\mu}(\mu, \rho, s) - \mu", dpi = 180,
		clims = (-μlim, μlim)
	)

	ρlim = maximum(abs.(dρ))
	
	g2fig = contourf(
		μspace, ρspace, dρ';
		c = :coolwarm, xlabel = L"f", ylabel = L"\rho", title = L"G_{\rho}(f, \rho, s) - \rho", dpi = 180,
		clims = (-ρlim, ρlim)
	)
	
	plot(g1fig, g2fig; size = (1000, 400), margins = 5Plots.mm)
end

# ╔═╡ c28509b7-df3e-44ba-b130-7626590f2f01
xbasin = ybasin = range(0.01, 0.99; length = 351)

# ╔═╡ 6dcc1295-1fe4-40a5-b691-fc392cf75c41
rvalues = range(1e-5, 0.5; length = 351)

# ╔═╡ 381fc3b1-55b3-48dc-85cd-1a4efab77d97
md"

## Dynamical System agent

``r`` $(@bind rfield Slider(
	0:0.001:1, 
	show_value = true, default = 1/8
))

"

# ╔═╡ dc3b28b0-def2-4c24-8104-51b9081d51e6
model = VerticalModel(m, 0.01, rfield, K)

# ╔═╡ 98904f31-4676-40a1-ab7c-da9c0b6911e0
md"
``\mu_0`` $(@bind μ₀ Slider(
	0:0.01:0.99, show_value = true, default = 0.5
))

``\rho_0`` $(@bind ρ₀ Slider(
	0:0.01:0.99, show_value = true, default = 0.3
))
"

# ╔═╡ ae7a7c1d-c3b6-4dc8-aeab-ef2af11ac25d
begin	

	X = Array{Float64}(undef, 2, model.K, 2) # (comp, soc), k, (f, ρ)
	X[1, 1, :] = X[2, 1, :] = [μ₀, ρ₀] 

	for k ∈ 2:model.K
		X[1, k, :] = G̃(X[1, k - 1, :], model; integer = false)
		X[2, k, :] = Gₚ(X[2, k - 1, :], m, rfield, μ₀, k)
	end

	vecfig = agentvectorfieldplot(model.m, model.r; L= 10, rescale = 4, integer = false)
	

	# Competitive
	plot!(vecfig, X[1, :, 1], X[1, :, 2], c = :red, label = "Agent", marker = :o, markersize = 2) 
	scatter!(vecfig, [X[1, 1, 1]], [X[1, 1, 2]], c = :red, label = nothing)

	if false
		# Social planner
		plot!(vecfig, X[2, :, 1], X[2, :, 2], c = :green, label = "Planner", marker = :o, markersize = 2)
		scatter!(vecfig, [X[2, 1, 1]], [X[2, 1, 2]], c = :green, label = nothing)
	end

	vecfig
end

# ╔═╡ 48b46e7a-5fc1-4b58-a9e8-73c2c532a3f9
md"
## Behaviour on the $\rho = 0$ axis.

``r`` $(@bind r0 Slider(
	range(0, 1 / ℯ + 0.3, length = 501), 
	show_value = true, default = 1 / 8.41
))

``\mu_0`` $(@bind initμ Slider(
	range(0, 1, length = 1001), 
	show_value = true, default = 0
))

"

# ╔═╡ 5a6baa64-9695-46cb-b59b-f289b0872ff7
begin
	function plotevolution(xs, fn::Function, args...; kwargs...)
		figure = plot(); plotevolution!(figure, xs, fn, args...; kwargs...)
		return figure
	end
	function plotevolution!(fig, xs, fn::Function, args...; kwargs...)
		xlims = extrema(xs)
				
		plot!(
			fig, xs, fn, 
			aspect_ratio = 1, xlims = xlims, ylims = xlims, 
			args...; kwargs...
		)
		
		plot!(fig, xs, x -> x, linestyle = :dash, c = :black, label = nothing)
		
		try 
			roots = find_zeros(x -> fn(x) - x, xlims...)
			colors = [abs(∂(fn, x₀)) < 1 ? :green : :red for x₀ ∈ roots]
			scatter!(fig, roots, fn.(roots), label = nothing, c = colors)
		catch end
		
	end

	# Cobweb plot
	"""
	Overlay the cobweb dynamics on function plot
	"""
	function cobwebplot(xs, fn::Function, x₀, T::Int, args...; kwargs...)
		figure = plot()
		cobwebplot!(figure, xs, fn, x₀, T, args...; kwargs...)
		return figure
	end
	function cobwebplot!(figure, xs, fn::Function, x₀, T::Int, args...; kwargs...)
		
		plotevolution!(figure, xs, fn; label = get(kwargs, :label, L"f_a(x)"), xlabel = L"x_t", ylabel = L"x_{t+1}", kwargs...)
		
		x, y = x₀, fn(x₀)
		scatter!([x], [y], label = nothing, color = :black)
		for t in 2:T
			plot!(figure, [x, y], [y, y], color = :black, label = nothing)
			plot!(figure, [y, y], [y, fn(y)], color = :black, label = nothing)
			x, y = y, fn(y)
		end
	end
	
end

# ╔═╡ 23fbfae9-f653-4352-a238-55a27beabe99
begin
	sint(μ; r) = floor((log(r) - log(-log(μ))) / log(μ))
	
	function g(μ; r)
		if μ ≈ 0 return 0. end
		if μ ≈ 1. return 1. end

		μ′ = -r/log(μ)
		
		return μ′ > 1 - 1e-3 ? 1. : μ′
	end

	
	gint(μ; r) = μ^sint(μ; r)
	
	g′(μ; r) = r / (μ * log(μ)^2)

	μ̄₀(r) = find_zeros(x -> x * log(x) + r, (0, 1))
end

# ╔═╡ e93c0097-0291-4c4e-bb3d-949d3e3ea3c6
begin
	εmargin = 0.01
	
	denseunit = range(0, 1; length = 501)
	steady_states = find_zeros(x -> x * log(x) + r0, (0, 1))
	rl = 1 / (ℯ - 0.1)

	gstable(μ) = g(μ; r = r0)
	gunstable(μ) = g(μ; r = rl)
	lowlabel = latexstring("\$ \\pi \\approx $(round(inv(rl), digits = 2)) \\kappa \$")
	highlabel = latexstring("\$ \\pi \\approx $(round(inv(r0), digits = 2)) \\kappa \$")
	
	evolfig = cobwebplot(
		denseunit, gstable, initμ, K;
		label = highlabel,
		c = :darkgreen, legend = :bottomright,
	)

	stable = [
		μ̄ ≤ r0 ? :darkblue : :darkred for μ̄ ∈ steady_states
	]


	cobwebplot!(
		evolfig, denseunit, gunstable, initμ, K;
		label = lowlabel,
		c = :darkorange,
		xlabel = L"\mu_k", ylabel = L"\mu_{k + 1}",
		margins = 5Plots.mm, xlims = (-εmargin, 1 + εmargin),  ylims = (-εmargin, 1 + εmargin),

	)

		scatter!(
		evolfig, steady_states, steady_states; 
		label = nothing, c = stable, markersize = 5
	)



	evolfig

end

# ╔═╡ 0b6a227a-d904-4519-b61d-683b92ecccfa
if SAVE savefig(evolfig, joinpath("../docs/plots", "one-dim-cobweb")) end

# ╔═╡ fcac7287-a664-4c44-aba4-cb3ba607a558
begin
	divr = ℯ^(-1) + 0.009
	convr = ℯ^(-1) - 0.009
	
	K₀ = 25
	X₀ = Matrix{Float64}(undef, (K₀, 2))
	X₀[1, :] .= 0.2
	
	for k ∈ 2:K₀
		X₀[k, 1] = g(X₀[k - 1, 1], r = convr)
		X₀[k, 2] = g(X₀[k - 1, 2], r = divr)
	end

	onedimtraj = plot(
		xlabel = latexstring("Layer, \$k\$"),
		ylabel = L"\mu_k", margins = 5Plots.mm, legend = :topleft,
		ylims = (0, 1 + εmargin),
		xticks = 0:(K₀ - 1), legendtitle = L"\kappa / \pi"
	)

	layers = 0:(K₀ - 1)
	plot!(onedimtraj, layers, k -> X₀[k + 1, 1]; marker = :o, c = :darkgreen, label = L"e^{-1} - \varepsilon")
	plot!(onedimtraj, layers, k -> X₀[k + 1, 2]; marker = :o, c = :darkorange, label = L"e^{-1} + \varepsilon")

	if SAVE savefig(onedimtraj, joinpath("../docs/plots", "one-dim-trajectory")) end

	onedimtraj

end

# ╔═╡ 59aff6b5-e80f-4aa4-99e0-21bc2b111434
begin
	rbifurcation = range(0, 1, length = 501)
	colors = [:darkblue, :darkred]
	orbitfig = vline(
		[1 / ℯ]; 
		linestyle = :dash, c = :black,
		aspect_ratio = 1 / 2, ylims = (-.01, 1.01), xlims = (-.01, 0.5),
		label = nothing, xlabel = L"\kappa / \pi", ylabel = "Equilibria",
		dpi = 180, legend = :topright, margins = 5Plots.mm
	)

	labels = ["stable" "unstable"]


	annotate!(
		orbitfig, 
		1/ℯ + 1e-2, 0.11, 
		text(
			latexstring("\$ e^{-1} \$"), :black, 18, :left
		)
	)
	
	for r ∈ rbifurcation

		
		μ̄ᵣ = μ̄₀(r)
		
		stabilitycolors = [colors[μ < r ? 1 : 2] for μ ∈ μ̄ᵣ]

		for (idx, μ) ∈ enumerate(μ̄ᵣ)
			label = r == 0.2 ? labels[idx] : nothing
			
			scatter!(
				orbitfig,
				[r], [μ];
				label = label, 
				c = stabilitycolors[idx],
				markersize = 1.5
			)
		end
	end

	orbitfig
end

# ╔═╡ 2fc1f58d-bf22-497b-ab3f-4c28d95b5bf9
if SAVE savefig(orbitfig, joinpath("../docs/plots", "one-dim-bif")) end

# ╔═╡ fc984dae-fa9f-43e1-abf9-9637563941fe
md"
### Basin of attraction

"

# ╔═╡ 64ae50a7-00af-4c09-b363-7ac686096a50
rbasin = 1 / 8

# ╔═╡ 779daf41-9752-4653-b3f0-9fe7de089396
function G̃!(dx, x, p, t)
	m, r = p
	model = VerticalModel(m, 0.01, r, K)
	dx .= G̃(x, model)
end

# ╔═╡ 25c861f8-9103-4e4f-a76f-b1816d842e31
ds = DiscreteDynamicalSystem(G̃!, [0.01, 0.01], [100, rbasin])

# ╔═╡ a20c0d6e-0e25-40d5-90bf-e156c4721b85
begin
	mapper = AttractorsViaRecurrences(ds, (xbasin, ybasin))
	basins, attractors = basins_of_attraction(mapper; show_progress = false)

	resattractors = Dict([k => v[1][1] for (k, v) in attractors])
	""
end

# ╔═╡ b5255c42-96a8-4c46-877c-18babf4e414a
begin
	foc(μ, ρ) = μ * (ψ₀(1 / ρ) - ψ₀(1 - μ + μ / ρ)) - rbasin
	interior_boundaries = find_zeros(μ -> μ * log(μ) + rbasin, (0.01, 0.99))

	boundary(μ) = find_zero(ρ -> foc(μ, ρ), μ)
end

# ╔═╡ d7ffac82-c55c-4e04-ac3b-efe8ad1d9fc6
begin
	conv = (b -> b > 0 ? 1 - resattractors[b] : 0).(basins)

	
	basinfig = contourf(
		xbasin, ybasin, conv'; 
		c = :viridis,
		aspect_ratio = 1, linewidth = 0, dpi = 180, 
		xticks = 0:0.1:1, yticks = 0:0.1:1,
		xlims = extrema(xbasin), ylims = extrema(ybasin),
		xlabel = L"\mu", ylabel = L"\rho", 
		levels = 0:0.01:1, clims = (0, 1)
	)

	annotate!(basinfig,
		1.2, 0.5, 
		text(
			latexstring("Attractors' resilience, \$1 - \\mu\$"), 
			:black, 18, rotation = -90
		)
	)

	plot!(
		basinfig,
		range(interior_boundaries..., length = 501), boundary;
		c = :white, linewidth = 3, label = nothing
	)
	
	basinfig
end

# ╔═╡ e901e035-1c0f-48ad-aae9-aee946027a04
if SAVE savefig(basinfig, joinpath("../docs/plots", "basin_small")) end

# ╔═╡ 91b18abe-f770-4613-9c84-8ebe9041ec98
begin

	rschoice = 1 / 8

	S̃space = ((t) -> agentoptimum(t[1], t[2]; m, r = rschoice)).(
		product(unit, unit)
	)

	extremaS = (0, ceil(maximum(S̃space)))
	
	agentsfig = contourf(
		unit, unit, S̃space',
		c = :viridis, linewidth = 0,
		clims = extremaS,
		levels = range(extremaS..., step = 0.1),
		dpi = 180,
		xlims = extrema(unit), ylims = extrema(unit),
		xlabel = L"\mu", ylabel = L"\rho",
		xticks = 0:0.1:1, yticks = 0:0.1:1,
		aspect_ratio = 1
	)

	contour!(
		agentsfig, unit, unit, 
		(f, ρ) -> agentoptimum(f, ρ; m, r = rschoice),
		linewidth = 3,
		levels = [1], c = :white
	)

	annotate!(
		agentsfig, 
		0.45, 0.75, text(L"\tilde{s} = 1", :white, FONT_SIZE)
	)

	annotate!(agentsfig,
		1.2, 0.5, 
		text(
			L"\tilde{s} \ (\mu, \rho)", 
			:black, 18, rotation = -90
		)
	)

	
	agentsfig
end

# ╔═╡ c2057c0c-b8f6-4741-bbf9-a22267e1c645
if SAVE savefig(agentsfig, joinpath("../docs/plots", "agents")) end

# ╔═╡ 777e44ed-5795-43c9-b90d-4bd1019f46ea
md"### Bifurcation"

# ╔═╡ 2d0b4c04-f29c-4fbd-ba29-c9b1a312c0f2
begin
	Nbif = 501
	paramspace = range(1 / 2ℯ, 1.2 / ℯ; length = Nbif)
	
	dsbif = DiscreteDynamicalSystem(G̃!, [0.4, 0.1], [100, 0.2])
	output = orbitdiagram(dsbif, [1, 2], 2, paramspace; n = Norbit, Ttr = Ttr)
end

# ╔═╡ 7cab141b-40cd-4906-91f2-9e33c19b9985
begin

	μbars = []
	ρbars = []

	for out ∈ output

		push!(μbars, out[1][1])
		push!(ρbars, out[1][2])
	end

	first_zero = findfirst(μ -> μ > 1 - 1e-3 , μbars)
	tipping_point = !isnothing(first_zero) ? first_zero : length(μbars)

	fig = plot(
		paramspace[1:(tipping_point - 1)], μbars[1:(tipping_point - 1)];
		c = :darkblue,
		ylabel = L"\mu",
		yguidefontcolor = :darkblue,
		xlabel = L"\kappa / \pi",
		legend = false, rightmargin = 25Plots.mm,
		alpha = 0.8, dpi = 300,
		ylims = (-0.01, 1.01),
		margins = 5Plots.mm
	)

	plot!(
		fig, paramspace[tipping_point:end], μbars[tipping_point:end],
		c = :darkblue
	)

	rtip = paramspace[tipping_point] 

	vline!(fig, [rtip], linestyle = :dash)

	plot!(
		twinx(fig), 
		paramspace[1:(tipping_point - 1)], ρbars[1:(tipping_point - 1)];
		ylabel = L"\rho",
		legend = false, c = :darkred,
		alpha = 0.8,
		yguidefontcolor = :darkred)

	plot!(
		fig, paramspace[tipping_point:end], ρbars[tipping_point:end],
		c = :darkred
	)

	
end

# ╔═╡ 33a48a25-d7b7-4947-bf31-0403e5b781f4
if SAVE savefig(fig, joinpath("../docs/plots", "bifurcation")) end

# ╔═╡ 782b9bd9-f78a-41f3-ad89-8fd222cfcac2
md"## Stable manifold for small $\rho$"

# ╔═╡ 5021c671-3515-4c94-96e1-797f287ce2f5
md"## Properties of the stable interior

``\rho`` $(@bind ρ̄ Slider(
	range(1e-3, 0.3; length = 101), 
	show_value = true, default = 1e-2
))
"

# ╔═╡ fc09995e-e8ca-47dc-a290-ef29390107ee
ρ̄ratio = (1 - ρ̄) / ρ̄

# ╔═╡ 06943ab9-186e-48bf-8f1e-f803a1df833d
foc(μ) = μ * (ψ₀(ρ̄ratio) - ψ₀(μ * ρ̄ratio) + (1 - μ) / ρ̄ratio)

# ╔═╡ b53ae8ed-f0c1-4e51-8ae5-44f755afaf6a
approximatefoc(μ) = ρ̄ * (μ^2 - μ/2 + 1/2) - μ * log(μ) 

# ╔═╡ 8125c155-1f02-4d25-9411-81e9cbb3a34e
let
	focfig = plot(unit, foc; xlabel = L"\mu", label = "Exact", title = "LHS first order condition")	
	
	plot!(focfig, unit, approximatefoc, label = "Linear approximation")
end

# ╔═╡ 8dddad86-fe0d-4421-b625-cfd1c7e8385c
function ρmanifold(μ; r)

	x = - (r + μ * log(μ)) / (μ^2 - μ/2 + 1/2)

	return x / (1 + x)
	
end

# ╔═╡ d85e6321-7afa-49aa-9cce-cd87f7db7960
begin
	rmanifoldplot = 1/3
	
	
	plot(0:0.01:1, μ -> ρmanifold(μ; r = rmanifoldplot), ylims = (0, 1))

end

# ╔═╡ ee127449-b7b3-4538-8db9-8c8bbafa4f38
md"
## Welfare analysis

``r`` $(@bind rwelfare Slider(
	range(0, 1 / ℯ + 0.3, length = 501), 
	show_value = true, default = 1 / (ℯ + 0.01)
))

``q`` $(@bind q Slider(
	range(0, 1; length = 101), 
	show_value = true, default = 0.01
))

"

# ╔═╡ 30239771-7353-4352-aa93-45fea4e31b2b
let
	N = 301
	modelwelfare = VerticalModel(100, 0.01, rwelfare, 100)
	
	contourf(
		range(0.01, 0.99; length = N), range(0.01, 0.99; length = N),
		(μ, ρ) -> W̃(q, μ, ρ; model = modelwelfare),
		linewidth = 0, clims = (0, 1)
	)
	
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
DynamicalSystems = "61744808-ddfa-5f27-97ff-6e42cc95d634"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
IterTools = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[compat]
Colors = "~0.12.8"
Combinatorics = "~1.0.2"
Distributions = "~0.25.64"
DynamicalSystems = "~2.3.0"
ForwardDiff = "~0.10.30"
Interpolations = "~0.13.6"
IterTools = "~1.4.0"
LaTeXStrings = "~1.3.0"
NLsolve = "~4.5.1"
Optim = "~1.7.0"
Plots = "~1.31.1"
PlutoUI = "~0.7.39"
Roots = "~2.0.1"
SpecialFunctions = "~2.1.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "Static"]
git-tree-sha1 = "6ccb71b40b04ad69152f1f83d5925de13911417e"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.19"

[[ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "7d255eb1d2e409335835dc8624c35d97453011eb"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.14"

[[ArrayInterfaceOffsetArrays]]
deps = ["ArrayInterface", "OffsetArrays", "Static"]
git-tree-sha1 = "c49f6bad95a30defff7c637731f00934c7289c50"
uuid = "015c0d05-e682-4f19-8f0a-679ce4c54826"
version = "0.1.6"

[[ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceStaticArraysCore", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "efb000a9f643f018d5154e56814e338b5746c560"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.4"

[[ArrayInterfaceStaticArraysCore]]
deps = ["Adapt", "ArrayInterfaceCore", "LinearAlgebra", "StaticArraysCore"]
git-tree-sha1 = "a1e2cf6ced6505cbad2490532388683f1e88c3ed"
uuid = "dd5226c6-a4d4-4bc7-8575-46859f9c95b9"
version = "0.1.0"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "eaee37f76339077f86679787a71990c4e465477f"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.4"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "b1a532a582dd18b34543366322d390e1560d40a9"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.23"

[[CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "32abd86e3c2025db5172aa182b982debed519834"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.1"

[[CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "2dd813e5f2f7eec2d1268c57cf2373d3ee91fcea"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.1"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[ChaosTools]]
deps = ["Clustering", "Combinatorics", "DSP", "DelayEmbeddings", "Distances", "Distributions", "DynamicalSystemsBase", "Entropies", "ForwardDiff", "IntervalRootFinding", "LinearAlgebra", "LombScargle", "Neighborhood", "ProgressMeter", "Random", "Roots", "SpecialFunctions", "StaticArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "85e8370059f8fcbd99d437d86bbac472cbab9a98"
uuid = "608a59af-f2a3-5ad4-90b4-758bdf3122a7"
version = "2.9.0"

[[CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "5522c338564580adf5d58d91e43a55db0fa5fb39"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.10"

[[Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "75479b7df4167267d75294d14b58244695beb2ac"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.2"

[[ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonSolve]]
git-tree-sha1 = "332a332c97c7071600984b3c31d9067e1a4e6e25"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.1"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "9be8be1d8a6f44b96482c8af52238ea7987da3e3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.45.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "59d00b3139a9de4eb961057eabb65ac6522be954"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.0"

[[Contour]]
git-tree-sha1 = "a599cfb8b1909b0f97c5e1b923ab92e1c0406076"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.1"

[[CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "32d125af0fb8ec3f8935896122c5e345709909e5"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.0"

[[Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "3fb5d9183b38fdee997151f723da42fb83d1c6f2"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.6"

[[DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelayEmbeddings]]
deps = ["Distances", "Distributions", "LinearAlgebra", "Neighborhood", "Random", "StaticArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "c067ea7738b5d1d2680831c1fd155e33e91d5dc1"
uuid = "5732040d-69e3-5649-938a-b6b4f237613f"
version = "2.2.0"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[DiffEqBase]]
deps = ["ArrayInterfaceCore", "ChainRulesCore", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "Printf", "RecursiveArrayTools", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "f7a479aac5f3917b8472ac5f1b77d6f296fe58f1"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.92.3"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "28d605d9a0ac17118fe2c5e9ce0fbb76c3ceb120"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.11.0"

[[Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "d530092b57aef8b96b27694e51c575b09c7f0b2e"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.64"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[DynamicalSystems]]
deps = ["ChaosTools", "DelayEmbeddings", "DynamicalSystemsBase", "Entropies", "RecurrenceAnalysis", "Reexport", "Scratch"]
git-tree-sha1 = "5f6afe8c904f66e2f74ad23f032d526e690a3f55"
uuid = "61744808-ddfa-5f27-97ff-6e42cc95d634"
version = "2.3.0"

[[DynamicalSystemsBase]]
deps = ["DelayEmbeddings", "ForwardDiff", "LinearAlgebra", "SciMLBase", "SimpleDiffEq", "SparseArrays", "StaticArrays", "Statistics"]
git-tree-sha1 = "4846191e79594f281f72173db80e0d683126f3d7"
uuid = "6e36e845-645a-534a-86f2-f5d4aa5a06b4"
version = "2.7.3"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Entropies]]
deps = ["DelayEmbeddings", "Distances", "LinearAlgebra", "Neighborhood", "SparseArrays", "SpecialFunctions", "StaticArrays", "Statistics", "Wavelets"]
git-tree-sha1 = "6da8090e8e49cebf107647b988b536492b8deb74"
uuid = "ed8fcbec-b94c-44b6-89df-898894ad9591"
version = "1.1.2"

[[ErrorfreeArithmetic]]
git-tree-sha1 = "d6863c556f1142a061532e79f611aa46be201686"
uuid = "90fa49ef-747e-5e6f-a989-263ba693cf1a"
version = "0.5.2"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[FastBroadcast]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "LinearAlgebra", "Polyester", "Static", "StrideArraysCore"]
git-tree-sha1 = "21cdeff41e5a1822c2acd7fc7934c5f450588e00"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.1"

[[FastRounding]]
deps = ["ErrorfreeArithmetic", "LinearAlgebra"]
git-tree-sha1 = "6344aa18f654196be82e62816935225b3b9abe44"
uuid = "fa42c844-2597-5d31-933b-ebd51ab2693f"
version = "0.3.1"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "ee13c773ce60d9e95a6c6ea134f25605dce2eda3"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.13.0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "2f18915445b248731ec5db4e4a17e451020bf21e"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.30"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "4078d3557ab15dd9fe6a0cf6f65e3d4937e98427"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c98aea696662d09e215ef7cda5296024a9646c75"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.4"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "3a233eeeb2ca45842fe100e0413936834215abf5"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.4+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "db5c7e27c0d46fd824d470a3c32a4fc6c935fa96"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.1"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "b7b88a4716ac33fe31d6556c02fc60017594343c"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.8"

[[HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "cb7099a0109939f16a4d3b572ba8396b1f6c7c31"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.10"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "a8671d5c9670a62cb36b7d44c376bdb09181aa26"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.3"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[IntervalArithmetic]]
deps = ["CRlibm", "FastRounding", "LinearAlgebra", "Markdown", "Random", "RecipesBase", "RoundingEmulator", "SetRounding", "StaticArrays"]
git-tree-sha1 = "421f305e970dd1d2c8339c93b7674fd3a698ed06"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.20.6"

[[IntervalRootFinding]]
deps = ["ForwardDiff", "IntervalArithmetic", "LinearAlgebra", "Polynomials", "Reexport", "StaticArrays"]
git-tree-sha1 = "b6969692c800cc5b90608fbd3be83189edc5e446"
uuid = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
version = "0.5.10"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "f3c7f871d642d244e7a27e3fb81e8441e13230d8"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.8.0"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "46a39b9c58749eefb5f2dc1178cb8fab5332b1ab"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.15"

[[LayoutPointers]]
deps = ["ArrayInterface", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "b67e749fb35530979839e7b4b606a97105fe4f1c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.10"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LombScargle]]
deps = ["FFTW", "LinearAlgebra", "Measurements", "Random", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d64a0ce7539181136a85fd8fe4f42626387f0f26"
uuid = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
version = "1.0.3"

[[LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SIMDTypes", "SLEEFPirates", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "7bf979d315193570cc2b79b4d2eb4595d68b9352"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.119"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "891d3b4e8f8415f53108b4918d0183e61e18015b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.0"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "dd8b9e6d7be9731fdaecc813acc5c3083496a251"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.7.2"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "4e675d6e9ec02061800d6cfb695812becbd03cdf"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.4"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "0e353ed734b1747fc20cd4cba0edd9ac027eff6a"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.11"

[[Neighborhood]]
deps = ["Distances", "NearestNeighbors", "Random", "Test"]
git-tree-sha1 = "1159fcaf3b72923cf623b2748f238a5115ed2623"
uuid = "645ca80c-8b79-4109-87ea-e1f58159d116"
version = "0.2.3"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[NonlinearSolve]]
deps = ["ArrayInterfaceCore", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "8a00c7b9418270f1fa57da319d11febbe5f92101"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.20"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "1ea784113a6aa054c5ebd95945fa5e52c2f378e7"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.7"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9a36165cf84cff35851809a40a928e1103702013"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.16+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "7a28efc8e34d5df89fc87343318b0a8add2c4021"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ca433b9e2f5ca3a0ce6702a032fce95a3b6e1e48"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.14"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9888e59493658e476d3073f1ce24348bdc086660"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.0"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "93e82cebd5b25eb33068570e3f63a86be16955be"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.31.1"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "97bbf8dc886d67ff0dd1f56cfc0ee18b7bb7f8ce"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.6.13"

[[PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "4cd738fca4d826bef1a87cbe43196b34fa205e6d"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.6"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "a1f7f4e41404bed760213ca01d7f384319f717a5"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.25"

[[PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[RecurrenceAnalysis]]
deps = ["DelayEmbeddings", "DelimitedFiles", "Distances", "Graphs", "LinearAlgebra", "Random", "SparseArrays", "StaticArrays", "Statistics", "UnicodePlots"]
git-tree-sha1 = "6ef236edbec1e6e48ef64d2bb41e571267213f3a"
uuid = "639c3291-70d9-5ea2-8c5b-839eba1ee399"
version = "1.8.1"

[[RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "ZygoteRules"]
git-tree-sha1 = "7ddd4f1ac52f9cc1b784212785f86a75602a7e4b"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.31.0"

[[RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "3ee71214057e29a8466f5d70cfe745236aa1d9d7"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.11"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[Roots]]
deps = ["CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "30e3981751855e2340e9b524ab58c1ec85c36f33"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.1"

[[RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "dd4195d308df24f33fb10dde7c22103ba88887fa"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.1"

[[SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "7ee0e13ac7cd77f2c0e93bff8c40c45f05c77a5a"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.33"

[[SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "RecipesBase", "RecursiveArrayTools", "StaticArraysCore", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "3243a883fa422a0a5cfe2d3b6ea6287fc396018f"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.42.2"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SetRounding]]
git-tree-sha1 = "d7a25e439d07a17b7cdf97eecee504c50fedf5f6"
uuid = "3cc68bcd-71a2-5612-b932-767ffbe40ab0"
version = "0.2.1"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleDiffEq]]
deps = ["DiffEqBase", "MuladdMacro", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "f8711f4e31bc8c10e59fd698bef155ab9278a50a"
uuid = "05bca326-078c-5bf0-a5bf-ce7c7982d7fd"
version = "1.5.1"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "a9e798cae4867e3a41cae2dd9eb60c047f1212db"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.6"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "46638763d3a25ad7818a15d441e0c3446a10742d"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.7.5"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "9f8a5dc5944dc7fbbe6eb4180660935653b0a9d9"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.0"

[[StaticArraysCore]]
git-tree-sha1 = "66fe9eb253f910fe8cf161953880cfdaef01cdf0"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.0.1"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2c11d7290036fe7aac9038ff312d3b3a2a5bf89e"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.4.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "48598584bacbebf7d30e20880438ed1d24b7c7d6"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.18"

[[StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "ac730bd978bf35f9fe45daa0bd1f51e493e97eb4"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.3.15"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "ec47fb6069c57f1cee2f67541bf8f23415146de7"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.11"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "f8629df51cab659d70d2e5618a430b4d3f37f2c3"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.0"

[[TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "0a4d8838dc28b4bcfaa3a20efb8d63975ad6781d"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.8.0"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "caf797b6fccbc0d080c44b4cb2319faf78c9d058"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.12"

[[Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[UnicodePlots]]
deps = ["Crayons", "Dates", "SparseArrays", "StatsBase"]
git-tree-sha1 = "3cb994143aba28cfe66615702505b2d294cebd3e"
uuid = "b8865327-cd53-5732-bb35-84acbb429228"
version = "2.5.1"

[[Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "39e55018bccc5a858217db32aa3d9e7decbefd0c"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.40"

[[Wavelets]]
deps = ["DSP", "FFTW", "LinearAlgebra", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "52e87ecea56e02e0672c7f3c9fd9ca03915d7e1b"
uuid = "29a6e085-ba6d-5f35-a997-948ac2efa89a"
version = "0.9.4"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─d3add18c-4845-4683-a82a-cbbe16f32b6f
# ╠═8b32f51d-e1f7-489d-95af-5e25909d709d
# ╠═963896ed-ae03-4b68-bec2-6b9917093454
# ╠═1a39488c-43c8-47e3-a760-913962cb3703
# ╠═c87845ef-2a91-41c8-933a-b6e139739927
# ╟─edf0c6f6-db46-11ec-19d0-f3a901cbdf5e
# ╟─82d52cc7-c6d4-45e5-8eeb-e81fda6317f8
# ╠═551ba67b-5317-4dd0-9b9d-126ca8c13eb0
# ╠═077b426a-ac05-4b6d-bcf1-0e379a82ba76
# ╠═7bb7dd41-0768-476f-b815-536607ecbdf0
# ╟─7253f132-6b20-43aa-8df8-fd8ef15cbcf5
# ╠═a0ff3a4e-c189-4459-b610-1a7109ffbcbc
# ╟─ab5e7adf-8acc-4fff-92f5-c5d66d23f7f7
# ╟─7c2dafe1-13ce-4c57-8c83-ef0d4b235302
# ╠═e8397cfc-407a-451d-b571-d95dd66c8c33
# ╟─15f7b2fb-f98c-464d-9055-59bb796adb7b
# ╠═1e41fda9-6e2e-4402-b685-72d49e6caf9a
# ╠═5c1f36ac-2d44-4011-af8c-7634473ffc54
# ╟─b67eb0ba-6db7-40ba-a0b2-21151c4eb748
# ╠═49e1c995-831b-472e-b68c-b56a21d1f308
# ╠═e518339b-9781-4f37-895f-19343315c79c
# ╟─5228005b-19f6-4cab-bd31-31aedb4fcf6c
# ╠═c28509b7-df3e-44ba-b130-7626590f2f01
# ╠═6dcc1295-1fe4-40a5-b691-fc392cf75c41
# ╟─381fc3b1-55b3-48dc-85cd-1a4efab77d97
# ╠═9f11ee9e-e89f-41c4-8ef4-d91c7d4e8db3
# ╠═dc3b28b0-def2-4c24-8104-51b9081d51e6
# ╟─98904f31-4676-40a1-ab7c-da9c0b6911e0
# ╟─ae7a7c1d-c3b6-4dc8-aeab-ef2af11ac25d
# ╠═48b46e7a-5fc1-4b58-a9e8-73c2c532a3f9
# ╟─5a6baa64-9695-46cb-b59b-f289b0872ff7
# ╟─23fbfae9-f653-4352-a238-55a27beabe99
# ╠═e93c0097-0291-4c4e-bb3d-949d3e3ea3c6
# ╠═0b6a227a-d904-4519-b61d-683b92ecccfa
# ╠═fcac7287-a664-4c44-aba4-cb3ba607a558
# ╠═59aff6b5-e80f-4aa4-99e0-21bc2b111434
# ╠═2fc1f58d-bf22-497b-ab3f-4c28d95b5bf9
# ╟─fc984dae-fa9f-43e1-abf9-9637563941fe
# ╠═64ae50a7-00af-4c09-b363-7ac686096a50
# ╠═779daf41-9752-4653-b3f0-9fe7de089396
# ╟─25c861f8-9103-4e4f-a76f-b1816d842e31
# ╠═a20c0d6e-0e25-40d5-90bf-e156c4721b85
# ╟─b5255c42-96a8-4c46-877c-18babf4e414a
# ╟─d7ffac82-c55c-4e04-ac3b-efe8ad1d9fc6
# ╠═e901e035-1c0f-48ad-aae9-aee946027a04
# ╠═91b18abe-f770-4613-9c84-8ebe9041ec98
# ╠═c2057c0c-b8f6-4741-bbf9-a22267e1c645
# ╟─777e44ed-5795-43c9-b90d-4bd1019f46ea
# ╠═2d0b4c04-f29c-4fbd-ba29-c9b1a312c0f2
# ╟─7cab141b-40cd-4906-91f2-9e33c19b9985
# ╠═33a48a25-d7b7-4947-bf31-0403e5b781f4
# ╟─782b9bd9-f78a-41f3-ad89-8fd222cfcac2
# ╟─5021c671-3515-4c94-96e1-797f287ce2f5
# ╠═fc09995e-e8ca-47dc-a290-ef29390107ee
# ╠═06943ab9-186e-48bf-8f1e-f803a1df833d
# ╠═b53ae8ed-f0c1-4e51-8ae5-44f755afaf6a
# ╟─8125c155-1f02-4d25-9411-81e9cbb3a34e
# ╠═8dddad86-fe0d-4421-b625-cfd1c7e8385c
# ╠═d85e6321-7afa-49aa-9cce-cd87f7db7960
# ╟─ee127449-b7b3-4538-8db9-8c8bbafa4f38
# ╠═30239771-7353-4352-aa93-45fea4e31b2b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
