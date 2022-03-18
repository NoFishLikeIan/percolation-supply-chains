function goodlayout(m::Vector{Int64})    

    g = length(m)

    xs = vcat(
        (mᵢ -> range(0, 1; length = mᵢ)).(m)...
    )

    ys = repeat(range(0, 1; length = g), inner = g)

    xs, ys

end

function plotgraph(goods::Goods, functional::Vector{Bool}, x::Vector{Vector{Int64}})

    m = length.(goods)
    n = sum(m)

    firms = vcat(goods...)

    g = SimpleDiGraph(n)
    for i ∈ firms, j ∈ x[i]
        Graphs.add_edge!(g, j, i)
    end

    # Colors

    goodcolors = distinguishable_colors(
        length(goods), [colorant"blue", colorant"red"]
    )

    colors = [
        functional[i] ? goodcolors[findfirst(g -> i ∈ g, goods)] : RGB(1, 1, 1)
        for i ∈ firms
    ]

    # Layout    
    xs, ys = goodlayout(m)

    gplot(g, xs, ys; nodelabel = firms, nodefillc = colors)

end