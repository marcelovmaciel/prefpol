

function boxplot_bootstraps(d;
        variants::Vector{String} = ["zero","random","mice"],
        ncols::Int = 3,
        nrows::Int = 2,
        figsize_px::Tuple = (1200,800))

    # helper: map variant => 1,2,3 …
    variant_idx(v) = findfirst(==(v), variants)

    measures = sort(collect(keys(d)))                # stable order
    fig      = Figure(resolution = figsize_px)

    for (i, m) in enumerate(measures)
        row = fld(i-1, ncols) + 1
        col = (i-1) % ncols + 1

        ax = Axis(fig[row, col];
                title  = string(m),
                xlabel = "imputation variant",
                ylabel = "value")

        xs = Int[]
        ys = Float32[]

        for v in variants
            vals = d[m][Symbol(v)]
            append!(xs, fill(variant_idx(v), length(vals)))
            append!(ys, Float32.(vals))
        end

        boxplot!(ax, xs, ys)
        ax.xticks = (1:length(variants), variants)
    end

    return fig
end




function boxplot_overlay(d;
    variants::Vector{String} = ["zero","random","mice"],
    palette ::NTuple = (:dodgerblue3, :darkorange2, :seagreen4),
    boxwidth::Real = 0.18,
    figsize::Tuple = (1200, 500))

    measures = sort(collect(keys(d)))          # stable order on x‑axis
nm, nv   = length(measures), length(variants)

# x‑offsets so the three boxes sit side‑by‑side at each integer tick
offsets = LinRange(-boxwidth, boxwidth, nv)

fig = Figure(resolution = figsize)
ax  = Axis(fig[1,1]; xlabel = "measure", ylabel = "value")

for (j, v) in enumerate(variants)          # one colour per variant
    xs, ys = Float32[], Float32[]

    for (i, m) in enumerate(measures)      # gather all y's for this variant
        vals = d[m][Symbol(v)]
        append!(xs, fill(i + offsets[j], length(vals)))
        append!(ys, Float32.(vals))
    end

    boxplot!(ax, xs, ys;
             color = palette[j],
             width = boxwidth * 0.9,
             label = v)
end

ax.xticks = (1:nm, string.(measures))
axislegend(ax; position = :rt)
return fig
end



function boxplot_cases_by_measure(d;
    variants::Vector{String} = ["zero","random","mice"],
    palette ::Vector = Makie.wong_colors(),
    boxwidth::Real = 0.18,
    figsize::Tuple = (1000, 500), n_alternatives =2, n_bootstrap = 5000)

measures = sort(collect(keys(d)))         # legend order
nv, nm   = length(variants), length(measures)

# offsets so nm boxes sit side‑by‑side at each variant tick
offsets = LinRange(-boxwidth, boxwidth, nm)

fig = Figure(resolution = figsize)
Label(fig[0, 1:2];
          text     = "B = $n_bootstrap • m = $n_alternatives",
          fontsize = 16,
          halign   = :center)
ax  = Axis(fig[1, 1]; xlabel = "imputation variant", ylabel = "value")

for (j, m) in enumerate(measures)         # iterate measures → colour
xs, ys = Float32[], Float32[]

for (i, v) in enumerate(variants)     # gather variant values
vals = d[m][Symbol(v)]
append!(xs, fill(i + offsets[j], length(vals)))
append!(ys, Float32.(vals))
end

col = palette[(j - 1) % length(palette) + 1]
labels = Dict(
    :calc_reversal_HHI             => "HHI",
    :calc_total_reversal_component => "R",
    :fast_reversal_geometric       => "RHHI"
)

    #= :weighted_psi_symmetric        => L"\Psi_{w}",
        :psi_we => L"\Psi_{we}",
        :psi_wb => L"\Psi_{wb}" =#
boxplot!(ax, xs, ys; color = col, 
width = boxwidth * 0.9, label = get(labels, m, string(m)))
end

ax.xticks = (1:nv, variants)
Legend(fig[1, 2], ax)   

return fig
end

 









#fig = plot_divergence_heatmap_cairo(divergence_by_religion, :D10)


function plot_divergence_heatmap_cairo(df::DataFrame, demo::Symbol;
    colormap = :viridis)

labels  = string.(df[!, demo])                               # row / col names
mat     = Matrix{Float64}(df[:, Not([demo, :proportion])])   # numeric part
n       = length(labels)

fig = Figure(resolution = (600, 600))
ax  = Axis(fig[1, 1]; aspect = DataAspect())

hm = heatmap!(ax, mat;                         # ← only the matrix
colormap    = colormap,
colorrange  = (minimum(mat), maximum(mat)),
interpolate = false)
# another option for colorrange would be colorrange =  (0.,1.)

# ticks exactly at the cell centres 1,2,…,n
ax.xticks = (1:n, labels)
ax.yticks = (1:n, labels)
ax.xticklabelrotation = π/4

ax.title  = "Pair-wise Divergence — $(demo)"
ax.xlabel = "Consensus from Group"
ax.ylabel = "Profiles from Group"

Colorbar(fig[1, 2], hm; label = "Divergence", width = 15)

return fig
end







"""
    boxplot_C_D_by_variant(stats::Dict;
                           variants = ["zero","random","mice"],
                           palette  = Makie.wong_colors(),
                           boxwidth = 0.18,
                           figsize  = (700, 400),
                           title    = "")

Produce a side-by-side box-plot of the **C** and **D** vectors stored in
`stats` (as produced by `bootstrap_group_metrics`).  
Each variant (:zero, :random, :mice) appears on the x-axis; the two measures
get different colours and horizontal offsets.

Returns a `Figure`.
"""
function boxplot_C_D_by_variant(stats::Dict;
                                variants::Vector{String}      = ["zero","random","mice"],
                                palette ::Vector              = Makie.wong_colors(),
                                boxwidth::Real                = 0.18,
                                figsize::Tuple{<:Integer,<:Integer} = (700, 400),
                                title::AbstractString         = "")

    measures = [:C, :D]                     # fixed two measures
    nm       = length(measures)
    nv       = length(variants)
    offsets  = LinRange(-boxwidth, boxwidth, nm)

    fig = Figure(resolution = figsize)
    ax  = Axis(fig[1, 1]; xlabel = "imputation variant", ylabel = "value")
    title != "" && (ax.title = title)

    for (j, m) in enumerate(measures)                    # colour loop
        xs, ys = Float32[], Float32[]
        for (i, v) in enumerate(variants)                # variant loop
            vals = stats[Symbol(v)][m]                   # Vector{Float64}
            append!(xs, fill(i + offsets[j], length(vals)))
            append!(ys, Float32.(vals))
        end
        col = palette[(j - 1) % length(palette) + 1]
        lab = m == :C ? "C (coherence)" : "D (divergence)"
        boxplot!(ax, xs, ys; color = col,
                 width = boxwidth*0.9, label = lab)
    end

    ax.xticks = (1:nv, variants)
    Legend(fig[1, 2], ax)
    return fig
end








"""
    combined_C_D_boxplots(race_stats, religion_stats, sex_stats;
                          variants  = ["zero","random","mice"],
                          titles    = ["Race (D12a)", "Religion (D10)", "Sex (D02)"],
                          palette   = Makie.wong_colors(),
                          boxwidth  = 0.18,
                          figsize   = (2100, 500))

Create a 1 × 3 panel of box-plots (one per grouping variable).  
Each panel shows the bootstrap distributions of **C** and **D** across the
imputation variants.  A single legend is placed on the right-hand side.
"""
function combined_C_D_boxplots(race_stats::Dict,
                               religion_stats::Dict,
                               sex_stats::Dict,
                               ideology_stats::Dict;
                               variants  = ["zero","random","mice"],
                               titles    = ["Race (D12a)", "Religion (D10)",
                                            "Sex (D02)",  "Ideology (Q19)"],
                               palette   = Makie.wong_colors(),
                               boxwidth  = 0.18,
                               figsize   = (1400, 900))

    stats_vec = [race_stats, religion_stats, sex_stats, ideology_stats]

    nv        = length(variants)
    measures  = [:C, :D]
    measure_labels = ["C (coherence)", "D (divergence)"]
    offsets   = LinRange(-boxwidth, boxwidth, length(measures))

    fig = Figure(resolution = figsize)

    # ── create 4 axes in a 2×2 grid ──────────────────────────────────
    axes = [Axis(fig[r, c];
                 title  = titles[(r-1)*2 + c],
                 xlabel = "imputation variant",
                 ylabel = "value")
            for r in 1:2, c in 1:2] |> vec

    # ── draw box-plots and capture two handles for legend ────────────
    legend_handles = BoxPlot[]; legend_labels = String[]

    for (p, stats) in enumerate(stats_vec)          # panel index 1–4
        ax = axes[p]

        for (j, m) in enumerate(measures)           # C / D
            xs = Float32[]; ys = Float32[]
            for (i, v) in enumerate(variants)       # zero / random / mice
                vals = stats[Symbol(v)][m]
                append!(xs, fill(i + offsets[j], length(vals)))
                append!(ys, Float32.(vals))
            end
            col = palette[(j-1) % length(palette) + 1]
            bp  = boxplot!(ax, xs, ys;
                           color = col,
                           width = boxwidth*0.9)

            if p == 1                               # take handles once
                push!(legend_handles, bp)
                push!(legend_labels, measure_labels[j])
            end
        end
        ax.xticks = (1:nv, variants)
    end

    # legend occupies the third column spanning both rows
    Legend(fig[1:2, 3], legend_handles, legend_labels)

    return fig
end





"""
    boxplot_alt_by_variant(measures_over_m;
                           variants   = ["zero","random","mice"],
                           palette    = Makie.wong_colors(),
                           boxwidth   = 0.18,
                           figsize    = (1000, 900))

`measures_over_m` is the Dict produced by `compute_measures_over_alternatives`,
mapping   m → Dict(measure → Dict(variant → Vector)).

Creates a 3×1 stacked figure (one row per variant),  
x-axis = number of alternatives, y-axis = value, coloured boxes = measures.
"""
function boxplot_alt_by_variant(measures_over_m::AbstractDict;
                                variants   ::Vector{String} = ["zero","random","mice"],
                                palette    ::Vector         = Makie.wong_colors(),
                                boxwidth   ::Real           = 0.18,
                                figsize    ::Tuple          = (1000, 900))

    # ——— gather keys --------------------------------------------------------
    ms        = sort(collect(keys(measures_over_m)))                    # e.g. [2,3,4,5]
    first_d   = first(values(measures_over_m))
    measures  = sort(collect(keys(first_d)))                            # e.g. [:C,:D,:G]

    nm, nv    = length(measures), length(variants)
    offsets   = LinRange(-boxwidth, boxwidth, nm)

    # legend labels (add more if you need)
    labels = Dict(
        :calc_reversal_HHI             => "HHI",
        :calc_total_reversal_component => "R",
        :fast_reversal_geometric       => "RHHI",
        :weighted_psi_symmetric        => L"\Psi_{w}",
        :psi_we => L"\Psi_{we}",
        :psi_wb => L"\Psi_{wb}"
    )

    # ——— figure & layout ----------------------------------------------------
    fig = Figure(resolution = figsize)
rowgap!(fig.layout, 18)          # <<<<<<  changed
colgap!(fig.layout, 8)           # <<<<<<  changed
fig[1, 2] = GridLayout()    
colsize!(fig.layout, 2, 180)     # legend column ≈180 px

    Label(fig[0, 1:2];
          text     = "Number of alternatives m = $(first(ms)) … $(last(ms))   •   bootstrap by variant",
          fontsize = 18,
          halign   = :center)

    axes = [Axis(fig[i, 1];
                 title  = variants[i],
                 xlabel = "number of alternatives",
                 ylabel = "value") for i in 1:nv]

    legend_handles = BoxPlot[]; legend_labels = AbstractString[]

    # ——— populate axes ------------------------------------------------------
    for (row, var) in enumerate(variants)
        ax = axes[row]

        for (j, meas) in enumerate(measures)
            xs = Float32[]; ys = Float32[]

            for m in ms
                vals = measures_over_m[m][meas][Symbol(var)]
                append!(xs, fill(Float32(m) + offsets[j], length(vals)))
                append!(ys, Float32.(vals))
            end

            col = palette[(j - 1) % length(palette) + 1]
            bp  = boxplot!(ax, xs, ys; color = col, width = boxwidth*0.9,
                                   label = get(labels, meas, string(meas)))

            if row == 1                      # collect legend once
                push!(legend_handles, bp)
                push!(legend_labels, get(labels, meas, string(meas)))
            end
        end
        ax.xticks = (ms, string.(ms))
    end

    Legend(fig[1:nv, 2], legend_handles, legend_labels)
    return fig
end

function describe_candidate_set(candidates::Vector{String})
    pretty_names = [join([uppercasefirst(lowercase(w)) for w in split(name, "_")], " ")
                    for name in candidates]
    return "Candidates: " * join(pretty_names, ", ")
end

"""
    lines_alt_by_variant(measures_over_m;
                         variants   = ["zero","random","mice"],
                         palette    = Makie.wong_colors(),
                         figsize    = (1000, 900))

`measures_over_m` is the Dict returned by
`compute_measures_over_alternatives`, i.e.

    m → measure → variant → Vector{Float64}

For every *variant* (row) and every *measure* (colour) it draws:
  • a line through the median of the bootstrap distribution  
  • a translucent band between the 25 % and 75 % quantiles.
"""
#= function lines_alt_by_variant(measures_over_m::AbstractDict;
                              variants   = ["zero","random","mice"],
                              palette    = Makie.wong_colors(),
                              figsize    = (1000, 900), candidate_label::String = "")

    ms        = sort(collect(keys(measures_over_m)))                # alt counts
    measures  = sort(collect(keys(first(values(measures_over_m))))) # :C, :D, …
    nm, nv    = length(measures), length(variants)

    # legend labels
    mlabels = Dict(
        :calc_reversal_HHI             => "HHI",
        :calc_total_reversal_component => "R",
        :fast_reversal_geometric       => "RHHI",
        :weighted_psi_symmetric        => L"\Psi_{w}",
    )

    # ───── figure & layout ──────────────────────────────────────────
    fig = Figure(resolution = figsize)
    rowgap!(fig.layout, 18)
    colgap!(fig.layout, 4)
    fig[1, 2] = GridLayout()                 # create legend column
    colsize!(fig.layout, 2, 100)

    first_m      = first(ms)
    first_meas   = first(measures)
    first_var    = first(variants)
    n_bootstrap  = length(measures_over_m[first_m][first_meas][Symbol(first_var)])

   #=  Label(fig[0, 1:2];
          text     = "Number of alternatives m = $(first(ms)) … $(last(ms))   •   $n_bootstrap bootstraps by variant",
          fontsize = 18,  halign   = :left,
      padding  = (0, 0, 0, 10))  =#
titlegrid = GridLayout(tellwidth =  true )

fig[0, 1] = titlegrid

Label(titlegrid[1, 1];
      text     = "Number of alternatives m = $(first(ms)) … $(last(ms))   •   $n_bootstrap bootstraps by variant",
      fontsize = 18,
      halign   = :left)

Label(titlegrid[2, 1];
      text     = candidate_label,
      fontsize = 14,
      halign   = :left)
#fig[0, 1] = titlegrid   

    axes = [Axis(fig[i, 1];
                 title        = variants[i],
                 xlabel       = "number of alternatives",
                 ylabel       = "value",
                 yticks  = (0:0.1:1, string.(0:0.1:1)) , 
                 xticks       = (ms, string.(ms)))
                 for i in 1:nv]

    legend_handles = Vector{Lines}()
    legend_labels  = AbstractString[]

    # ───── plot per row / per measure ───────────────────────────────
    for (row, var) in enumerate(variants)
        ax = axes[row]

        for (j, meas) in enumerate(measures)
            col  = palette[(j - 1) % length(palette) + 1]
            meds = Float64[]                       # medians per m
            q25s = Float64[]                       # 25 % quantile
            q75s = Float64[]                       # 75 % quantile

            for m in ms
                vals = measures_over_m[m][meas][Symbol(var)]
                push!(meds, median(vals))
                push!(q25s, quantile(vals, 0.25))
                push!(q75s, quantile(vals, 0.75))
            end

            # band (IQR)
            band!(ax, ms, q25s, q75s;
                  color = (col, 0.25),             # add alpha
                  linewidth = 0)

            # median line
            ln = lines!(ax, ms, meds;
                        color = col, linewidth = 2,
                        label = get(mlabels, meas, string(meas)))

            # collect legend entry once
            if row == 1
                push!(legend_handles, ln)
                push!(legend_labels, get(mlabels, meas, string(meas)))
            end
        end
    end

    Legend(fig[1:nv, 2], legend_handles, legend_labels)
    resize_to_layout!(fig)  
    return fig
end
 =#


function lines_alt_by_variant(measures_over_m::AbstractDict;
                              variants   = ["zero","random","mice"],
                              palette    = Makie.wong_colors(),
                              figsize    = (1000, 900),
                              candidate_label::String = "", year)

    ms        = sort(collect(keys(measures_over_m)))                # alt counts
    measures  = sort(collect(keys(first(values(measures_over_m))))) # :C, :D, …
    nv        = length(variants)

    # legend labels
    mlabels = Dict(
        :calc_reversal_HHI             => "HHI",
        :calc_total_reversal_component => "R",
        :fast_reversal_geometric       => "RHHI",
    )
#=     :weighted_psi_symmetric        => L"\Psi_{w}",
        :psi_we => L"\Psi_{we}",
        :psi_wb => L"\Psi_{wb}"
 =#
    # ── figure & layout ─────────────────────────────────────────────
    fig = Figure(resolution = figsize)
    rowgap!(fig.layout, 18)
    colgap!(fig.layout, 4)
    fig[1, 2] = GridLayout()                 # legend column
    colsize!(fig.layout, 2, 100)

    # header
    first_m, last_m = first(ms), last(ms)
    first_var = first(variants)
    n_bootstrap = length(first(values(first(values(measures_over_m))))[Symbol(first_var)])

    titlegrid = GridLayout(tellwidth = true)
    fig[0, 1] = titlegrid
    Label(titlegrid[1, 1];
          text = "Year = $(year)   •   Number of alternatives = $first_m … $last_m   •   $n_bootstrap pseudo-profiles",
          fontsize = 18, halign = :left)
    Label(titlegrid[2, 1];
          text = candidate_label,
          fontsize = 14, halign = :left)

    # axes
    axes = [Axis(fig[i, 1];
                 title   = variants[i],
                 xlabel  = "number of alternatives",
                 ylabel  = "value",
                 yticks  = (0:0.1:1, string.(0:0.1:1)),
                 xticks  = (ms, string.(ms)))
            for i in 1:nv]

    legend_handles = Lines[]
    legend_labels  = AbstractString[]

    # ── plot ────────────────────────────────────────────────────────
    for (row, var) in enumerate(variants)
        ax = axes[row]

        for (j, meas) in enumerate(measures)
            col = palette[(j-1) % length(palette) + 1]

            meds = Float64[]          # medians per m
            q25s = Float64[]; q75s = Float64[]          # 25–75%
            p05s = Float64[]; p95s = Float64[]          # 5–95%  ### NEW

            for m in ms
                vals = measures_over_m[m][meas][Symbol(var)]
                push!(meds, median(vals))
                push!(q25s, quantile(vals, 0.25))
                push!(q75s, quantile(vals, 0.75))
                push!(p05s, quantile(vals, 0.05))      ### NEW
                push!(p95s, quantile(vals, 0.95))      ### NEW
            end

            # 90 % band  (5‒95)   ### NEW
            band!(ax, ms, p05s, p95s;
                  color = (col, 0.12), linewidth = 0)

            # inter-quartile band (25‒75)
            band!(ax, ms, q25s, q75s;
                  color = (col, 0.25), linewidth = 0)

            # median line
            ln = lines!(ax, ms, meds;
                        color = col, linewidth = 2,
                        label = get(mlabels, meas, string(meas)))

            # legend once
            if row == 1
                push!(legend_handles, ln)
                push!(legend_labels, get(mlabels, meas, string(meas)))
            end
        end
    end

    Legend(fig[1:nv, 2], legend_handles, legend_labels)
    resize_to_layout!(fig)
    return fig
end



nx = PyCall.pyimport("networkx")  


function _margin_graph_py(stats;
                          title="Margin graph",
                          digits=2,
                          figsize=(10, 8))
nx = PyCall.pyimport("networkx") 
    G = nx.DiGraph()
    edge_labels = Dict{Tuple{String,String},String}()

    # ── nodes & edges ─────────────────────────────────────────────────────
    for ((a,b), (μ,σ)) in stats
        abs(μ) < 1e-12 && continue
        G.add_node(string(a));  G.add_node(string(b))

        from,to = μ > 0 ? (a,b) : (b,a)
        lbl = string(round(abs(μ);digits=digits)," ± ",round(σ;digits=digits))
        G.add_edge(string(from), string(to))
        edge_labels[(string(from), string(to))] = lbl
    end

    # layout
    pos = nx.spring_layout(G, seed=42)

    # figure
    plt.figure(figsize=figsize)

    # nodes
    nx.draw_networkx_nodes(G, pos,
                           node_color="#DDDDDD",
                           edgecolors="black",
                           node_size=2500)
    nx.draw_networkx_labels(G, pos, font_size=10)

    # curved edges (all of them, single path)
    curved = G.edges()
    nx.draw_networkx_edges(G, pos,
                           edgelist=curved,
                           connectionstyle="arc3,rad=0.15",
                           arrowstyle="-|>",
                           arrows=true,
                           arrowsize=30,
                           edge_color="#555555",
                           width=1.,
                           min_source_margin=15,
                           min_target_margin=15)

    # edge labels (offset slightly perpendicular to edge)
    for (u,v) in curved
        x1,y1 = pos[u]
        x2,y2 = pos[v]
        xm, ym = (x1+x2)/2, (y1+y2)/2           # midpoint

        # perpendicular offset
        dx, dy = y2 - y1, -(x2 - x1)
        norm = sqrt(dx^2 + dy^2) + 1e-9
        dx, dy = dx/norm, dy/norm
        xm += 0.05*dx;   ym += 0.05*dy          # nudge

        plt.text(xm, ym, edge_labels[(u,v)],
                 fontsize=9, ha="center", va="center", color="black")
    end

    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    return plt.gcf()
end





#= 

function lines_group_measures_over_m(
    stats_by_m;                                   # Dict{m → …}  from your pipeline
    demographics::Vector{Symbol},
    variants::Vector{String} = ["zero","random","mice"],
    measures::Vector{Symbol} = [:C,:D,:G],
    m_values::Vector{Int}    = sort(collect(keys(stats_by_m))),
    palette                  = Makie.wong_colors(),
    maxcols::Int             = 3,
    n_yticks::Int            = 5,                # total y-ticks incl. min & max
    figsize                  = (300 * min(maxcols, length(demographics)),
                                300 * ceil(Int, length(demographics)/maxcols))
)

    # colour per measure, line-style per variant
    measure_cols   = Dict(measures[i] => palette[i] for i in eachindex(measures))
    variant_styles = Dict("zero"=>:solid, "random"=>:dash, "mice"=>:dot)

    # grid layout
    n_demo = length(demographics)
    ncol   = min(maxcols, n_demo)
    nrow   = ceil(Int, n_demo / ncol)
    fig    = Figure(resolution = figsize)
    rowgap!(fig.layout, 24); colgap!(fig.layout, 24)

    # legend collectors
    legend_handles = Any[]
    legend_labels  = String[]

    for (idx, demo) in enumerate(demographics)
        r, c = fldmod1(idx, ncol)

        ax = Axis(fig[r, c];
                  title  = string(demo),
                  xlabel = "number of alternatives",
                  ylabel = "value",
                  xticks = (m_values, string.(m_values)))

        # gather values to compute min/max later
        allvals = Float64[]

        for meas in measures, var in variants
            data_per_m = [stats_by_m[m][demo][Symbol(var)][meas] for m in m_values]
            append!(allvals, vcat(data_per_m...))

            meds = mean.(data_per_m)
            q25s = map(x -> quantile(x, 0.25), data_per_m)
            q75s = map(x -> quantile(x, 0.75), data_per_m)

            col = measure_cols[meas]
            sty = variant_styles[var]

            band!(ax, m_values, q25s, q75s; color = (col, 0.20), linewidth = 0)
            ln  = lines!(ax, m_values, meds; color = col, linestyle = sty, linewidth = 2)

            if idx == 1
                push!(legend_handles, ln)
                push!(legend_labels, "$(meas) • $var")
            end
        end

        # ---- add min & max to y-ticks ----------------------------------
        y_min, y_max = minimum(allvals), maximum(allvals)
        ticks = range(y_min, y_max; length = n_yticks) |> collect
        labels = string.(round.(ticks; digits = 3))
        ax.yticks[] = (ticks, labels)                 # overwrite ticks
    end

    Legend(fig[1:nrow, ncol+1], legend_handles, legend_labels; tellheight = false)
    resize_to_layout!(fig)
    return fig
end =#


#= function lines_group_measures_over_m(
    stats_by_m;
    year,
    demographics::Vector{Symbol},
    variants::Vector{String} = ["zero","random","mice"],
    measures::Vector{Symbol} = [:C,:D,:G],
    m_values::Vector{Int}    = sort(collect(keys(stats_by_m))),
    candidate_label::String  = "",
    palette                  = Makie.wong_colors(),
    maxcols::Int             = 3,
    n_yticks::Int            = 5,
    figsize                  = (300 * min(maxcols, length(demographics)),
                                300 * ceil(Int, length(demographics)/maxcols)),
)

    # ——— metadata for the title ————————————————————————————————
    first_m      = m_values[1]
    first_demo   = demographics[1]
    first_var    = variants[1]
    first_meas   = measures[1]
    n_boot       = length(stats_by_m[first_m][first_demo][Symbol(first_var)][first_meas])
    first_m, last_m = first(m_values), last(m_values)
    title_txt    = "Year = $(year) • $(n_boot) bootstraps • number of alternatives = $first_m … $last_m"
    subtitle_txt = isempty(candidate_label) ? "" : candidate_label

    # ——— colour & style dictionaries ————————————————————————————
    measure_cols   = Dict(measures[i] => palette[i] for i in eachindex(measures))
    variant_styles = Dict("zero"=>:solid, "random"=>:dash, "mice"=>:dot)

    # ——— layout bookkeeping (extra rows for title / subtitle) ————
    n_demo      = length(demographics)
    ncol        = min(maxcols, n_demo)
    nrow        = ceil(Int, n_demo / ncol)
    extra_rows  = 1 + (!isempty(subtitle_txt) ? 1 : 0)     # title [+ subtitle]
    row_shift   = extra_rows                               # axes start after these

    fig = Figure(resolution = figsize)
    rowgap!(fig.layout, 24);  colgap!(fig.layout, 24)

    # figure-level labels
    fig[1, 1:ncol] = Label(fig, title_txt;  fontsize = 20, halign = :left)
    if !isempty(subtitle_txt)
        fig[2, 1:ncol] = Label(fig, subtitle_txt; fontsize = 14, halign = :left)
    end

    # ——— legend collectors ————————————————————————————————
    legend_handles = Any[]
    legend_labels  = String[]

    # ——— main loop ————————————————————————————————————————————
    for (idx, demo) in enumerate(demographics)
        r, c = fldmod1(idx, ncol)                 # 1-based grid coords for demo
        ax = Axis(fig[r + row_shift, c];
                  title  = string(demo),
                  xlabel = "number of alternatives",
                  ylabel = "value",
                  xticks = (m_values, string.(m_values)))

        # collect data to fix y-ticks per panel
        allvals = Float64[]

        for meas in measures, var in variants
            data_per_m = [stats_by_m[m][demo][Symbol(var)][meas] for m in m_values]
            append!(allvals, vcat(data_per_m...))

            meds = mean.(data_per_m)
            q25s = map(x -> quantile(x, 0.25), data_per_m)
            q75s = map(x -> quantile(x, 0.75), data_per_m)

            col = measure_cols[meas]
            sty = variant_styles[var]

            band!(ax, m_values, q25s, q75s;
                  color = (col, 0.20), linewidth = 0)
            ln = lines!(ax, m_values, meds;
                         color = col, linestyle = sty, linewidth = 2)

            if idx == 1         # only populate legend once
                push!(legend_handles, ln)
                push!(legend_labels, "$(meas) • $var")
            end
        end

        # ——— tidy y-axis ————————————————————————————————
        y_min, y_max = extrema(allvals)
        ticks  = range(y_min, y_max; length = n_yticks) |> collect
        labels = string.(round.(ticks; digits = 3))
        ax.yticks[] = (ticks, labels)
    end

    # place legend at the right of the grid
    Legend(fig[row_shift+1 : row_shift+nrow, ncol+1],
           legend_handles, legend_labels; tellheight = false)

    resize_to_layout!(fig)
    return fig
end =#



function lines_group_measures_over_m(
    stats_by_m;                                    # Dict{m → …}
    year,                                          # Int or String
    demographics,                                  # Vector of Symbols/Strings
    variants         = ["zero","random","mice"],
    measures         = [:C,:D,:G],
    m_values         = nothing,                    # will default below
    candidate_label  = "",
    palette          = Makie.wong_colors(),
    maxcols          = 3,
    n_yticks         = 5,
    figsize          = nothing,                    # will default below
)

    # ── fill in defaults that depend on other kwargs ───────────────────
    m_values === nothing && (m_values = sort(collect(keys(stats_by_m))))
    n_demo  = length(demographics)

    if figsize === nothing
        figsize = (300 * min(maxcols, n_demo),
                   300 * ceil(Int, n_demo / maxcols))
    end

    # ── metadata for the title ─────────────────────────────────────────
    first_m, last_m = first(m_values), last(m_values)
    first_demo      = demographics[1]
    n_boot = length(stats_by_m[first_m][first_demo][Symbol(variants[1])][measures[1]])

    title_txt = "Year = $(year) • $(n_boot) pseudo-profiles • number of alternatives = $first_m … $last_m"
    subtitle_txt = isempty(candidate_label) ? "" : candidate_label

    # ── colour & style dictionaries ────────────────────────────────────
    measure_cols   = Dict(measures[i] => palette[i] for i in eachindex(measures))
    variant_styles = Dict("zero" => :solid, "random" => :dash, "mice" => :dot)

    # ── grid geometry ──────────────────────────────────────────────────
    ncol       = min(maxcols, n_demo)
    nrow       = ceil(Int, n_demo / ncol)
    extra_rows = 1 + (!isempty(subtitle_txt) ? 1 : 0)
    row_shift  = extra_rows

    # ── widen canvas if title lines are long ───────────────────────────
    est_title_px = 10 * max(length(title_txt), length(subtitle_txt))
    fig_width    = max(first(figsize), est_title_px + 60)
    fig_height   = last(figsize)
    fig          = Figure(resolution = (fig_width, fig_height))
    rowgap!(fig.layout, 24);  colgap!(fig.layout, 24)

    # ── figure-level labels ────────────────────────────────────────────
    fig[1, 1:ncol] = Label(fig, title_txt;  fontsize = 20, halign = :left)
    if !isempty(subtitle_txt)
        fig[2, 1:ncol] = Label(fig, subtitle_txt; fontsize = 14, halign = :left)
    end

    # ── legend collectors ──────────────────────────────────────────────
    legend_handles = Any[]
    legend_labels  = String[]

    # ── main loop over demographic panels ──────────────────────────────
    for (idx, demo) in enumerate(demographics)
        r, c = fldmod1(idx, ncol)
        ax = Axis(fig[r + row_shift, c];
                  title  = string(demo),
                  xlabel = "number of alternatives",
                  ylabel = "value",
                  xticks = (m_values, string.(m_values)))

        allvals = Float64[]

        for meas in measures, var in variants
            data_per_m = [stats_by_m[m][demo][Symbol(var)][meas] for m in m_values]
            append!(allvals, vcat(data_per_m...))

            meds = mean.(data_per_m)
            q25s = map(x -> quantile(x, 0.25), data_per_m)
            q75s = map(x -> quantile(x, 0.75), data_per_m)

            col = measure_cols[meas]
            sty = variant_styles[var]

            band!(ax, m_values, q25s, q75s; color = (col, 0.20), linewidth = 0)
            ln = lines!(ax, m_values, meds; color = col, linestyle = sty, linewidth = 2)

            if idx == 1
                push!(legend_handles, ln)
                push!(legend_labels, "$(meas) • $var")
            end
        end

        # nice y-ticks
        y_min, y_max = extrema(allvals)
        ticks  = range(y_min, y_max; length = n_yticks) |> collect
        labels = string.(round.(ticks; digits = 3))
        ax.yticks[] = (ticks, labels)
    end

    # legend at the right of the grid
    Legend(fig[row_shift + 1 : row_shift + nrow, ncol + 1],
           legend_handles, legend_labels; tellheight = false)

    resize_to_layout!(fig)
    return fig
end





const _PRETTY_MEASURE = Dict(
    :calc_reversal_HHI             => "HHI",
    :calc_total_reversal_component => "R",
    :fast_reversal_geometric       => "R HHI",
    :weighted_psi_symmetric        => "Ψ₍w₎",
    :psi_we                        => "Ψ₍we₎",
    :psi_wb                        => "Ψ₍wb₎",
    :C => "C",  :D => "D",  :G => "G"
)

function compare_global_measures_v3(
        res_dict;
        variants ::Vector{String} = ["mice"],
        palette  ::Vector         = Makie.wong_colors(),
        figsize  ::Tuple{<:Integer,<:Integer} = (1200, 850))

    @assert !isempty(variants) "variants list cannot be empty"

    scen_names = sort(collect(keys(res_dict)))
    first_res  = first(values(res_dict))

    ms        = sort(collect(keys(first_res.measures_over_m)))
    measures  = sort(collect(keys(first(values(first_res.measures_over_m)))))
    n_meas    = length(measures)

    year   = first_res.year
    n_boot = length(first(values(first(values(first_res.measures_over_m))))[Symbol(variants[1])])

    # ── grid geometry --------------------------------------------------------
    ncol = min(3, ceil(Int, sqrt(n_meas)))          # 4→2, 5–9→3
    nrow = ceil(Int, n_meas / ncol)

    fig = Figure(resolution = figsize)
    rowgap!(fig.layout, 18); colgap!(fig.layout, 22)

    fig[0, 1:ncol] = Label(fig,
        "Year $year   •   number of alternatives = $(first(ms))…$(last(ms))   •   $n_boot pseudo-profiles",
        fontsize = 20, halign = :left)

    # ── create axes ----------------------------------------------------------
    axes = Dict{Symbol,Axis}()
    for (idx, meas) in enumerate(measures)
        r = div(idx - 1, ncol) + 1
        c = mod(idx - 1, ncol) + 1
        axes[meas] = Axis(fig[r, c];
                          title  = get(_PRETTY_MEASURE, meas, string(meas)),
                          xlabel = "number of alternatives",
                          ylabel = "value",
                          xticks = (ms, string.(ms)))
    end

    # colour per scenario -----------------------------------------------------
    scen_col = Dict(name => palette[(i-1) % length(palette) + 1]
                    for (i, name) in enumerate(scen_names))

    legend_handles = Lines[]; legend_labels = String[]

    # ── drawing loop ---------------------------------------------------------
    for scen in scen_names
        colour = scen_col[scen]
        data   = res_dict[scen].measures_over_m
        first_line = nothing

        for meas in measures, var in variants
            ax   = axes[meas]
            meds = Float64[]; q25s=Float64[]; q75s=Float64[]
            p05s = Float64[]; p95s=Float64[]

            for m in ms
                vals = data[m][meas][Symbol(var)]
                push!(meds, median(vals))
                push!(q25s, quantile(vals, 0.25));  push!(q75s, quantile(vals, 0.75))
                push!(p05s, quantile(vals, 0.05));  push!(p95s, quantile(vals, 0.95))
            end

            # 90 % band (transparent)
            band!(ax, ms, p05s, p95s; color = (colour, 0.12), linewidth = 0)
            # IQR band
            band!(ax, ms, q25s, q75s; color = (colour, 0.25), linewidth = 0)

            ln = lines!(ax, ms, meds; color = colour, linewidth = 2)
            first_line === nothing && (first_line = ln)
        end

        push!(legend_handles, first_line)
        push!(legend_labels,
              describe_candidate_set(res_dict[scen].candidate_set))
    end

    # ── stacked legend at the bottom ----------------------------------------
    legend_row = nrow + 1
    Legend(fig[legend_row, 1:ncol],
           legend_handles, legend_labels;
           orientation = :vertical, framevisible = false, halign = :left)

    resize_to_layout!(fig)      # trims any unused whitespace
    return fig
end




#= 
function compare_demographic_across_scenarios(
        results_all::Dict,
        scenario_vec::Vector{Tuple{Int,String}};
        demographic,
        variant      ::String         = "mice",
        measures     ::Vector{Symbol} = [:C,:D,:G],
        palette      ::Vector         = Makie.wong_colors()[1:3],
        n_yticks     ::Int            = 5,
        figsize      ::Tuple{<:Integer,<:Integer} = (900, 320))

    @assert length(variant) > 0  "provide one imputation variant"
    demo_sym    = Symbol(demographic)
    variant_sym = Symbol(variant)
    n_panels    = length(scenario_vec)

    # —— meta from first scenario ————————————————————————————————
    first_year, first_scen = scenario_vec[1]
    stats_any  = results_all[first_year][first_scen].group_stats
    m_values   = sort(collect(keys(stats_any)))
    n_boot     = length(stats_any[first(m_values)][demo_sym][variant_sym][measures[1]])

    # —— figure & global title ————————————————————————————————
    fig = Figure(resolution = (figsize[1], figsize[2] * n_panels))
    rowgap!(fig.layout, 20)

    fig[0, 1] = Label(fig,
        "$demographic   •   m = $(first(m_values))…$(last(m_values))   •   B = $n_boot",
        fontsize = 20)

    # —— colour helpers ————————————————————————————————————————
    ΔL = Dict(:C => 0,    # base colour
              :D => +20,  # lighter
              :G => -15)  # darker
    RGB = Colors.RGB
    LCHab = Colors.LCHab
    meas_colour(base_rgb, meas) = begin
        base_lch = convert(LCHab, base_rgb)
        convert(RGB, LCHab(clamp(base_lch.l + ΔL[meas], 0, 100),
                           base_lch.c, base_lch.h))
    end

    legend_handles = Lines[]; legend_labels = ["C", "D", "G"]

    # —— panels loop ————————————————————————————————————————
    for (idx, (year, scen)) in enumerate(scenario_vec)
        stats       = results_all[year][scen].group_stats
        cand_label  = describe_candidate_set(results_all[year][scen].candidate_set)
        base_rgb    = convert(RGB, palette[(idx-1) % length(palette) + 1])

        ax = Axis(fig[idx, 1];
                  title  = "Year $year — $cand_label",
                  xlabel = "number of alternatives",
                  ylabel = "value",
                  xticks = (m_values, string.(m_values)))

        all_vals = Float64[]   # gather values for nice y-ticks

        for meas in measures
            col = meas_colour(base_rgb, meas)

            meds = Float64[]; q25 = Float64[]; q75 = Float64[]
            p05  = Float64[]; p95 = Float64[]

            for m in m_values
                vals = stats[m][demo_sym][variant_sym][meas]
                append!(all_vals, vals)

                push!(meds, median(vals))
                push!(q25, quantile(vals, 0.25));  push!(q75, quantile(vals, 0.75))
                push!(p05, quantile(vals, 0.05));  push!(p95, quantile(vals, 0.95))
            end

            band!(ax, m_values, p05, p95; color=(col, 0.12), linewidth=0)
            band!(ax, m_values, q25, q75; color=(col, 0.25), linewidth=0)
            ln = lines!(ax, m_values, meds; color=col, linewidth=2)

            idx == 1 && push!(legend_handles, ln)
        end

        # y-ticks based on all plotted values in this panel
        y_min, y_max = extrema(all_vals)
        ticks = range(y_min, y_max; length=n_yticks) |> collect
        ax.yticks[] = (ticks, string.(round.(ticks; digits=3)))
    end

    # —— one horizontal legend ————————————————————————————————
    Legend(fig[n_panels + 1, 1],
           legend_handles, legend_labels;
           orientation = :horizontal, framevisible = false)

    resize_to_layout!(fig)
    return fig
end =#
#= 

function compare_demographic_across_scenarios(
        results_all::Dict,
        scenario_vec::Vector{Tuple{Int,String}};
        demographic,
        variant      ::String         = "mice",
        measures     ::Vector{Symbol} = [:C,:D,:G],
        palette      ::Vector         = Makie.wong_colors()[1:3],
        n_yticks     ::Int            = 5,
        figsize      ::Tuple{<:Integer,<:Integer} = (900, 320))
    @assert length(variant) > 0  "provide one imputation variant"
    demo_sym    = Symbol(demographic)
    variant_sym = Symbol(variant)
    n_panels    = length(scenario_vec)
    # —— meta from first scenario ————————————————————————————————
    first_year, first_scen = scenario_vec[1]
    stats_any  = results_all[first_year][first_scen].group_stats
    m_values   = sort(collect(keys(stats_any)))
    n_boot     = length(stats_any[first(m_values)][demo_sym][variant_sym][measures[1]])
    # —— figure & global title ————————————————————————————————
    fig = Figure(resolution = (figsize[1], figsize[2] * n_panels))
    rowgap!(fig.layout, 20)
    fig[0, 1] = Label(fig,
        "$demographic   •   m = $(first(m_values))…$(last(m_values))   •   B = $n_boot",
        fontsize = 20)
    # —— colour helpers ————————————————————————————————————————
    ΔL = Dict(:C => 0,    # base colour
              :D => +20,  # lighter
              :G => -15)  # darker
    RGB = Colors.RGB
    LCHab = Colors.LCHab
    meas_colour(base_rgb, meas) = begin
        base_lch = convert(LCHab, base_rgb)
        convert(RGB, LCHab(clamp(base_lch.l + ΔL[meas], 0, 100),
                           base_lch.c, base_lch.h))
    end
    
    # —— assign consistent colors for measures ————————————————
    measure_colors = Dict()
    for (i, meas) in enumerate(measures)
        base_rgb = convert(RGB, palette[min(i, length(palette))])
        measure_colors[meas] = meas_colour(base_rgb, meas)
    end
    
    legend_handles = Lines[]; legend_labels = ["C", "D", "G"]
    # —— panels loop ————————————————————————————————————————
    for (idx, (year, scen)) in enumerate(scenario_vec)
        stats       = results_all[year][scen].group_stats
        cand_label  = describe_candidate_set(results_all[year][scen].candidate_set)
        ax = Axis(fig[idx, 1];
                  title  = "Year $year — $cand_label",
                  titlesize = 14,
                  titlegap = 8,
                  xlabel = "number of alternatives",
                  ylabel = "value",
                  xticks = (m_values, string.(m_values)))
        all_vals = Float64[]   # gather values for nice y-ticks
        for meas in measures
            col = measure_colors[meas]  # Use consistent color for this measure
            meds = Float64[]; q25 = Float64[]; q75 = Float64[]
            p05  = Float64[]; p95 = Float64[]
            for m in m_values
                vals = stats[m][demo_sym][variant_sym][meas]
                append!(all_vals, vals)
                push!(meds, median(vals))
                push!(q25, quantile(vals, 0.25));  push!(q75, quantile(vals, 0.75))
                push!(p05, quantile(vals, 0.05));  push!(p95, quantile(vals, 0.95))
            end
            band!(ax, m_values, p05, p95; color=(col, 0.12), linewidth=0)
            band!(ax, m_values, q25, q75; color=(col, 0.25), linewidth=0)
            ln = lines!(ax, m_values, meds; color=col, linewidth=2)
            idx == 1 && push!(legend_handles, ln)
        end
        # y-ticks based on all plotted values in this panel
        y_min, y_max = extrema(all_vals)
        ticks = range(y_min, y_max; length=n_yticks) |> collect
        ax.yticks[] = (ticks, string.(round.(ticks; digits=3)))
    end
    # —— one horizontal legend ————————————————————————————————
    Legend(fig[n_panels + 1, 1],
           legend_handles, legend_labels;
           orientation = :horizontal, framevisible = false)
    resize_to_layout!(fig)
    return fig
end =#

#= 
function compare_demographic_across_scenarios(
        results_all::Dict,
        scenario_vec::Vector{Tuple{Int,String}};
        demographic,
        variant      ::String         = "mice",
        measures     ::Vector{Symbol} = [:C,:D,:G],
        palette      ::Vector         = Makie.wong_colors()[1:3],
        n_yticks     ::Int            = 5,
        figsize      ::Tuple{<:Integer,<:Integer} = (900, 320))
    @assert length(variant) > 0  "provide one imputation variant"
    demo_sym    = Symbol(demographic)
    variant_sym = Symbol(variant)
    n_panels    = length(scenario_vec)
    # —— meta from first scenario ————————————————————————————————
    first_year, first_scen = scenario_vec[1]
    stats_any  = results_all[first_year][first_scen].group_stats
    m_values   = sort(collect(keys(stats_any)))
    n_boot     = length(stats_any[first(m_values)][demo_sym][variant_sym][measures[1]])
    # —— figure & global title ————————————————————————————————
    fig = Figure(resolution = (figsize[1], figsize[2] * n_panels))
    rowgap!(fig.layout, 20)
    fig[0, 1] = Label(fig,
        "$demographic   •   m = $(first(m_values))…$(last(m_values))   •   B = $n_boot",
        fontsize = 20)
    # —— colour helpers ————————————————————————————————————————
    ΔL = Dict(:C => 0,    # base colour
              :D => +20,  # lighter
              :G => -15)  # darker
    RGB = Colors.RGB
    LCHab = Colors.LCHab
    meas_colour(base_rgb, meas) = begin
        base_lch = convert(LCHab, base_rgb)
        convert(RGB, LCHab(clamp(base_lch.l + ΔL[meas], 0, 100),
                           base_lch.c, base_lch.h))
    end
    
    # —— assign consistent colors for measures ————————————————
    measure_colors = Dict()
    for (i, meas) in enumerate(measures)
        base_rgb = convert(RGB, palette[min(i, length(palette))])
        measure_colors[meas] = meas_colour(base_rgb, meas)
    end
    
    legend_handles = Lines[]; legend_labels = ["C", "D", "G"]
    # —— panels loop ————————————————————————————————————————
    for (idx, (year, scen)) in enumerate(scenario_vec)
        stats       = results_all[year][scen].group_stats
        cand_label  = describe_candidate_set(results_all[year][scen].candidate_set)
        # Wrap long titles by splitting at commas
        wrapped_title = replace(cand_label, ", " => ",\n")
        ax = Axis(fig[idx, 1];
                  title  = "Year $year — $wrapped_title",
                  titlesize = 12,
                  titlefont = :regular,  # Remove bold formatting
                  titlegap = 5,
                  xlabel = "number of alternatives",
                  ylabel = "value",
                  xticks = (m_values, string.(m_values)))
        all_vals = Float64[]   # gather values for nice y-ticks
        for meas in measures
            col = measure_colors[meas]  # Use consistent color for this measure
            meds = Float64[]; q25 = Float64[]; q75 = Float64[]
            p05  = Float64[]; p95 = Float64[]
            for m in m_values
                vals = stats[m][demo_sym][variant_sym][meas]
                append!(all_vals, vals)
                push!(meds, median(vals))
                push!(q25, quantile(vals, 0.25));  push!(q75, quantile(vals, 0.75))
                push!(p05, quantile(vals, 0.05));  push!(p95, quantile(vals, 0.95))
            end
            band!(ax, m_values, p05, p95; color=(col, 0.12), linewidth=0)
            band!(ax, m_values, q25, q75; color=(col, 0.25), linewidth=0)
            ln = lines!(ax, m_values, meds; color=col, linewidth=2)
            idx == 1 && push!(legend_handles, ln)
        end
        # y-ticks based on all plotted values in this panel
        y_min, y_max = extrema(all_vals)
        ticks = range(y_min, y_max; length=n_yticks) |> collect
        ax.yticks[] = (ticks, string.(round.(ticks; digits=3)))
    end
    # —— one horizontal legend ————————————————————————————————
    Legend(fig[n_panels + 1, 1],
           legend_handles, legend_labels;
           orientation = :horizontal, framevisible = false)
    resize_to_layout!(fig)
    return fig
end =#



#= function compare_demographic_across_scenarios(
        results_all::Dict,
        scenario_vec::Vector{Tuple{Int,String}};
        demographic,
        variant      ::String         = "mice",
        measures     ::Vector{Symbol} = [:C,:D,:G],
        palette      ::Vector         = Makie.wong_colors()[1:3],
        n_yticks     ::Int            = 5,
        figsize      ::Tuple{<:Integer,<:Integer} = (900, 320))
    @assert length(variant) > 0  "provide one imputation variant"
    demo_sym    = Symbol(demographic)
    variant_sym = Symbol(variant)
    n_panels    = length(scenario_vec)
    # —— meta from first scenario ————————————————————————————————
    first_year, first_scen = scenario_vec[1]
    stats_any  = results_all[first_year][first_scen].group_stats
    m_values   = sort(collect(keys(stats_any)))
    n_boot     = length(stats_any[first(m_values)][demo_sym][variant_sym][measures[1]])
    # —— figure & global title ————————————————————————————————
    fig = Figure(resolution = (figsize[1], figsize[2] * n_panels))
    rowgap!(fig.layout, 20)
    fig[0, 1:2] = Label(fig,
        "$demographic   •   m = $(first(m_values))…$(last(m_values))   •   B = $n_boot",
        fontsize = 20)
    # —— colour helpers ————————————————————————————————————————
    ΔL = Dict(:C => 0,    # base colour
              :D => +20,  # lighter
              :G => -15)  # darker
    RGB = Colors.RGB
    LCHab = Colors.LCHab
    meas_colour(base_rgb, meas) = begin
        base_lch = convert(LCHab, base_rgb)
        convert(RGB, LCHab(clamp(base_lch.l + ΔL[meas], 0, 100),
                           base_lch.c, base_lch.h))
    end
    
    # —— assign consistent colors for measures ————————————————
    measure_colors = Dict()
    for (i, meas) in enumerate(measures)
        base_rgb = convert(RGB, palette[min(i, length(palette))])
        measure_colors[meas] = meas_colour(base_rgb, meas)
    end
    
    legend_handles = Lines[]; legend_labels = ["C", "D", "G"]
    # —— panels loop ————————————————————————————————————————
    for (idx, (year, scen)) in enumerate(scenario_vec)
        stats       = results_all[year][scen].group_stats
        cand_label  = describe_candidate_set(results_all[year][scen].candidate_set)
        # Wrap long titles by splitting at commas
        wrapped_title = replace(cand_label, ", " => ",\n")
        ax = Axis(fig[1, idx];
                  title  = "Year $year — $wrapped_title",
                  titlesize = 12,
                  titlefont = :regular,  # Remove bold formatting
                  titlegap = 5,
                  xlabel = "number of alternatives",
                  ylabel = "value",
                  xticks = (m_values, string.(m_values)))
        all_vals = Float64[]   # gather values for nice y-ticks
        for meas in measures
            col = measure_colors[meas]  # Use consistent color for this measure
            meds = Float64[]; q25 = Float64[]; q75 = Float64[]
            p05  = Float64[]; p95 = Float64[]
            for m in m_values
                vals = stats[m][demo_sym][variant_sym][meas]
                append!(all_vals, vals)
                push!(meds, median(vals))
                push!(q25, quantile(vals, 0.25));  push!(q75, quantile(vals, 0.75))
                push!(p05, quantile(vals, 0.05));  push!(p95, quantile(vals, 0.95))
            end
            band!(ax, m_values, p05, p95; color=(col, 0.12), linewidth=0)
            band!(ax, m_values, q25, q75; color=(col, 0.25), linewidth=0)
            ln = lines!(ax, m_values, meds; color=col, linewidth=2)
            idx == 1 && push!(legend_handles, ln)
        end
        # y-ticks based on all plotted values in this panel
        y_min, y_max = extrema(all_vals)
        ticks = range(y_min, y_max; length=n_yticks) |> collect
        ax.yticks[] = (ticks, string.(round.(ticks; digits=3)))
    end
    # —— one horizontal legend ————————————————————————————————
    Legend(fig[n_panels + 1, 1:2],
           legend_handles, legend_labels;
           orientation = :horizontal, framevisible = false)
    resize_to_layout!(fig)
    return fig
end =#

#= function compare_demographic_across_scenarios(
        results_all::Dict,
        scenario_vec::Vector{Tuple{Int,String}};
        demographic::String,
        variant      ::String         = "mice",
        measures     ::Vector{Symbol} = [:C,:D,:G],
        palette      ::Vector         = Makie.wong_colors()[1:3],
        n_yticks     ::Int            = 5,
        base_width   ::Int            = 1000,      # a bit wider
        base_height  ::Int            = 360)       # a bit taller

    demo_sym    = Symbol(demographic)
    variant_sym = Symbol(variant)
    n_panels    = length(scenario_vec)

    # pull out m_values and B from first scenario
    y0, s0   = scenario_vec[1]
    stats0   = results_all[y0][s0].group_stats
    m_values = sort(collect(keys(stats0)))
    n_boot   = length(stats0[first(m_values)][demo_sym][variant_sym][:C])

    # helper to shade base colors per measure
    ΔL = Dict(:C=>0, :D=>+20, :G=>-15)
    function meas_colour(base_rgb, meas)
        lch = convert(Colors.LCHab, base_rgb)
        newL = clamp(lch.l + ΔL[meas], 0, 100)
        return convert(Colors.RGB, Colors.LCHab(newL, lch.c, lch.h))
    end

    # precompute one color per measure (so C,D,G stay consistent)
    base_rgbs = convert.(Colors.RGB, palette[1:length(measures)])
    measure_cols = Dict(measures[i] => meas_colour(base_rgbs[i], measures[i])
                        for i in eachindex(measures))

    # create figure: one row per panel + one for legend
    fig = Figure(resolution = (base_width, base_height * (n_panels + 1)))
    rowgap!(fig.layout, 30)      # more vertical space
    colgap!(fig.layout, 20)

    # global header
    fig[0, 1] = Label(fig,
        "$demographic   •   m=$(first(m_values))…$(last(m_values))   •   B=$n_boot",
        fontsize = 22,
        halign   = :left)

    legend_handles = Lines[]
    legend_labels  = String[]

    # loop through scenarios, stack vertically in col=1
    for (idx, (year, scen)) in enumerate(scenario_vec)
        stats    = results_all[year][scen].group_stats
        cand_lbl = describe_candidate_set(results_all[year][scen].candidate_set)
        long_title = "Year $year — $cand_lbl"

        # wrap at ~50 chars
        wrapped_title = join(TextWrap.wrap(long_title; width=50))

        ax = Axis(fig[idx, 1];
            title     = wrapped_title,
            titlesize = 14,
            titlegap  = 8,
            xlabel    = "number of alternatives",
            ylabel    = "value",
            xticks    = (m_values, string.(m_values)))

        allvals = Float64[]

        for meas in measures
            col  = measure_cols[meas]
            meds = Float64[]; q25 = Float64[]; q75 = Float64[]
            p05  = Float64[]; p95 = Float64[]

            for m in m_values
                vals = stats[m][demo_sym][variant_sym][meas]
                append!(allvals, vals)
                push!(meds, median(vals))
                push!(q25, quantile(vals,0.25)); push!(q75, quantile(vals,0.75))
                push!(p05, quantile(vals,0.05)); push!(p95, quantile(vals,0.95))
            end

            # draw 90% & IQR bands + median line
            band!(ax, m_values, p05, p95; color=(col,0.12), linewidth=0)
            band!(ax, m_values, q25, q75; color=(col,0.25), linewidth=0)
            ln = lines!(ax, m_values, meds; color=col, linewidth=2)

            # only grab legend handles on the first panel
            if idx == 1
                push!(legend_handles, ln)
                push!(legend_labels, string(meas))
            end
        end

        # tighten y-ticks to data range
        y_min, y_max = extrema(allvals)
        ticks = range(y_min, y_max; length=n_yticks) |> collect
        ax.yticks[] = (ticks, string.(round.(ticks; digits=3)))
    end

    # single legend in last row, col 1
    Legend(fig[n_panels+1, 1],
        legend_handles, legend_labels;
        orientation = :horizontal,
        framevisible = false,
        halign = :left)

    resize_to_layout!(fig)
    return fig
end =#

#= 
function compare_demographic_across_scenarios(
        results_all::Dict,
        scenario_vec::Vector{Tuple{Int,String}};
        demographic::String,
        variant      ::String         = "mice",
        measures     ::Vector{Symbol} = [:C,:D,:G],
        palette      ::Vector         = Makie.wong_colors()[1:3],
        n_yticks     ::Int            = 5,
        base_width   ::Int            = 400,       # width per panel
        base_height  ::Int            = 360)       # height for plots + header + legend

    # symbol conversions
    demo_sym    = Symbol(demographic)
    variant_sym = Symbol(variant)

    # how many panels
    n_panels = length(scenario_vec)

    # extract m and B from first scenario
    y0, s0   = scenario_vec[1]
    stats0   = results_all[y0][s0].group_stats
    m_values = sort(collect(keys(stats0)))
    n_boot   = length(stats0[first(m_values)][demo_sym][variant_sym][:C])

    # prepare a consistent colour for C/D/G
    ΔL = Dict(:C=>0, :D=>+20, :G=>-15)
    function shade(base_rgb::Colors.RGB, meas::Symbol)
        lch = convert(Colors.LCHab, base_rgb)
        newL = clamp(lch.l + ΔL[meas], 0, 100)
        convert(Colors.RGB, Colors.LCHab(newL, lch.c, lch.h))
    end

    # pick one base RGB per measure
    base_rgbs = convert.(Colors.RGB, palette[1:length(measures)])
    measure_cols = Dict(measures[i] => shade(base_rgbs[i], measures[i])
                        for i in eachindex(measures))

    # make figure: one row of panels + header row + legend row
    fig = Figure(
      resolution = (base_width * n_panels, base_height),
      layout = (3, n_panels)
    )
    rowgap!(fig.layout, 20)
    colgap!(fig.layout, 30)

    # global header spans all columns in row 1
    fig[1, 1:n_panels] = Label(fig,
        "$demographic • number of alternatives =$(first(m_values))…$(last(m_values)) • $n_boot bootstraps",
        fontsize = 22, halign = :center)

    # collect legend handles
    legend_handles = Lines[]
    legend_labels  = String[]

    # panel loop: place each at row 2, col i
    for (i, (year, scen)) in enumerate(scenario_vec)
        stats   = results_all[year][scen].group_stats
        cand_lbl = describe_candidate_set(results_all[year][scen].candidate_set)
        long_title = "Year $year — $cand_lbl"

        # wrap at ~50 chars
        wrapped = join(TextWrap.wrap(long_title; width=50))

        ax = Axis(fig[2, i];
            title     = wrapped,
            titlesize = 14,
            titlefont = "sans",
            titlegap  = 8,
            xlabel    = "number of alternatives",
            ylabel    = "value",
            xticks    = (m_values, string.(m_values))
        )

        allvals = Float64[]

        for meas in measures
            col  = measure_cols[meas]
            meds = Float64[]; q25 = Float64[]; q75 = Float64[]
            p05  = Float64[]; p95 = Float64[]

            for m in m_values
                vals = stats[m][demo_sym][variant_sym][meas]
                append!(allvals, vals)
                push!(meds, median(vals))
                push!(q25, quantile(vals, 0.25)); push!(q75, quantile(vals, 0.75))
                push!(p05, quantile(vals, 0.05)); push!(p95, quantile(vals, 0.95))
            end

            band!(ax, m_values, p05, p95; color=(col,0.12), linewidth=0)
            band!(ax, m_values, q25, q75; color=(col,0.25), linewidth=0)
            ln = lines!(ax, m_values, meds; color=col, linewidth=2)

            # only once, grab for the legend
            if i == 1
                push!(legend_handles, ln)
                push!(legend_labels, string(meas))
            end
        end

        # tighten y‐axis ticks
        y_min, y_max = extrema(allvals)
        ticks = range(y_min, y_max; length = n_yticks) |> collect
        ax.yticks[] = (ticks, string.(round.(ticks; digits=3)))
    end

    # legend under the panels, row 3
    fig[3, 1:n_panels] = Legend(fig,
        legend_handles, legend_labels;
        orientation  = :horizontal,
        framevisible = false,
        halign       = :center
    )

    resize_to_layout!(fig)
    return fig
end =#

function compare_demographic_across_scenariosy(
        results_all::Dict,
        scenario_vec::Vector{Tuple{Int,String}};
        demographic::String,
        variant      ::String         = "mice",
        measures     ::Vector{Symbol} = [:C,:D,:G],
        palette      ::Vector         = Makie.wong_colors()[1:3],
        n_yticks     ::Int            = 5,
        base_width   ::Int            = 400,       # width per panel
        base_height  ::Int            = 360)       # height for plots + header + legend

    # symbol conversions
    demo_sym    = Symbol(demographic)
    variant_sym = Symbol(variant)

    # how many panels
    n_panels = length(scenario_vec)

    # extract m and B from first scenario
    y0, s0   = scenario_vec[1]
    stats0   = results_all[y0][s0].group_stats
    m_values = sort(collect(keys(stats0)))
    n_boot   = length(stats0[first(m_values)][demo_sym][variant_sym][:C])

    # prepare a consistent colour for C/D/G
    ΔL = Dict(:C=>0, :D=>+20, :G=>-15)
    function shade(base_rgb::Colors.RGB, meas::Symbol)
        lch = convert(Colors.LCHab, base_rgb)
        newL = clamp(lch.l + ΔL[meas], 0, 100)
        convert(Colors.RGB, Colors.LCHab(newL, lch.c, lch.h))
    end

    # pick one base RGB per measure
    base_rgbs = convert.(Colors.RGB, palette[1:length(measures)])
    measure_cols = Dict(measures[i] => shade(base_rgbs[i], measures[i])
                        for i in eachindex(measures))

    # FIRST PASS: collect all values across all scenarios to determine global y-limits
    global_allvals = Float64[]
    
    for (year, scen) in scenario_vec
        stats = results_all[year][scen].group_stats
        for meas in measures
            for m in m_values
                vals = stats[m][demo_sym][variant_sym][meas]
                append!(global_allvals, vals)
            end
        end
    end
    
    # Calculate global y-limits
    global_y_min, global_y_max = extrema(global_allvals)
    global_ticks = range(global_y_min, global_y_max; length = n_yticks) |> collect

    # make figure: one row of panels + header row + legend row
    fig = Figure(
      resolution = (base_width * n_panels, base_height),
      layout = (3, n_panels)
    )
    rowgap!(fig.layout, 20)
    colgap!(fig.layout, 30)

    # global header spans all columns in row 1
    fig[1, 1:n_panels] = Label(fig,
        "$demographic • number of alternatives = $(first(m_values))…$(last(m_values)) • $n_boot pseudo-profiles",
        fontsize = 22, halign = :center)

    # collect legend handles
    legend_handles = Lines[]
    legend_labels  = String[]

    # panel loop: place each at row 2, col i
    for (i, (year, scen)) in enumerate(scenario_vec)
        stats   = results_all[year][scen].group_stats
        cand_lbl = describe_candidate_set(results_all[year][scen].candidate_set)
        long_title = "Year $year — $cand_lbl"

        # wrap at ~50 chars
        wrapped = join(TextWrap.wrap(long_title; width=50))

        ax = Axis(fig[2, i];
            title     = wrapped,
            titlesize = 14,
            titlefont = "sans",
            titlegap  = 8,
            xlabel    = "number of alternatives",
            ylabel    = "value",
            xticks    = (m_values, string.(m_values)),
            # Set shared y-limits and ticks
            limits    = (nothing, (global_y_min, global_y_max)),
            yticks    = (global_ticks, string.(round.(global_ticks; digits=3)))
        )

        for meas in measures
            col  = measure_cols[meas]
            meds = Float64[]; q25 = Float64[]; q75 = Float64[]
            p05  = Float64[]; p95 = Float64[]

            for m in m_values
                vals = stats[m][demo_sym][variant_sym][meas]
                push!(meds, median(vals))
                push!(q25, quantile(vals, 0.25)); push!(q75, quantile(vals, 0.75))
                push!(p05, quantile(vals, 0.05)); push!(p95, quantile(vals, 0.95))
            end

            band!(ax, m_values, p05, p95; color=(col,0.12), linewidth=0)
            band!(ax, m_values, q25, q75; color=(col,0.25), linewidth=0)
            ln = lines!(ax, m_values, meds; color=col, linewidth=2)

            # only once, grab for the legend
            if i == 1
                push!(legend_handles, ln)
                push!(legend_labels, string(meas))
            end
        end
    end

    # legend under the panels, row 3
    fig[3, 1:n_panels] = Legend(fig,
        legend_handles, legend_labels;
        orientation  = :horizontal,
        framevisible = false,
        halign       = :center
    )

    resize_to_layout!(fig)
    return fig
end
