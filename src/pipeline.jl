function build_bootstrapped_encoded_variants(scores_df::DataFrame;
                                             candidate_cols::Vector{String},
                                             demographics::Vector{String},
                                             n_bootstrap::Int = 500, ncandidates::Int = 5, 
                                             force_include::Vector{String}=String[])

    # 1. Compute candidate distributions & unique scores
    countmaps = build_candidate_score_distributions(scores_df, candidate_cols)
    # all_scores = extract_unique_scores(countmaps)
    countmaps2 = sanitize_countmaps(countmaps)
    nrespondents = nrow(scores_df)

    # 2. Compute "don't know her" metric and select top candidates
    dont_know_her = compute_dont_know_her(countmaps2, nrespondents)
    # most_known_candidates = get_most_known_candidates(dont_know_her, ncandidates)
    most_known_candidates = select_top_candidates(
        countmaps2, nrespondents;
        m = ncandidates,
        force_include = force_include)  
    println(most_known_candidates)
    # 3. Imputation variants using top known candidates
   
    top_imp = imputation_variants(scores_df, most_known_candidates, demographics;
                                  most_known_candidates = most_known_candidates)

    # 4. Encode imputation variants to code/profile mapping
    enc_top = encode_imputation_variants(top_imp, most_known_candidates, demographics)

    # 5. Bootstrap encoded variants
    bt_top = bootstrap_encoded_variants(enc_top;
                                        weights = scores_df.peso,
                                        B = n_bootstrap)

    return bt_top
end




function compute_bootstrap_measures(
        scores_df::DataFrame;
        candidate_cols::Vector{String} = CANDIDATOS_eseb2022,
        demographics::Vector{String}   = ["D02","D10","D12a","Q19_cat"],
        n_bootstrap::Int,
        n_alternatives::Int, 
        force_include::Vector{String}=String[]
    )

    bt_top = build_bootstrapped_encoded_variants(
        scores_df,
        candidate_cols  = candidate_cols,
        demographics    = demographics,
        n_bootstrap     = n_bootstrap,
        ncandidates     = n_alternatives, 
        force_include = force_include
    )

    return apply_all_measures_to_encoded_bts(bt_top)
end


function run_bootstrap_analysis(
        scores_df::DataFrame;
        n_bootstrap::Int,
        n_alternatives::Int,
        candidate_cols::Vector{String} = CANDIDATOS_eseb2022,
        demographics::Vector{String}   = ["D02","D10","D12a","Q19_cat"], force_include::Vector{String}=String[]
    )

    # get the measures
    measures_in_bts = compute_bootstrap_measures(
        scores_df;
        candidate_cols  = candidate_cols,
        demographics    = demographics,
        n_bootstrap     = n_bootstrap,
        n_alternatives  = n_alternatives, force_include = force_include)

    # plot
    fig = boxplot_cases_by_measure(
        measures_in_bts;
        n_bootstrap     = n_bootstrap,
        n_alternatives  = n_alternatives
    )

    # save
    fname = "exploratory/imgs/bootstrapTop$(n_alternatives)_cases_by_measure_B$(n_bootstrap).png"
    save(fname, fig; px_per_unit = 2)

    return fig
end



"""
    run_bootstrap_over_alternatives(
        scores_df::DataFrame,
        alt_list::Vector{Int};
        n_bootstrap::Int,
        candidate_cols::Vector{String}=CANDIDATOS_eseb2022,
        demographics::Vector{String}  = ["D02","D10","D12a","Q19_cat"]
    ) -> Dict{Int, Figure}

For each `m` in `alt_list`, calls `run_bootstrap_analysis(...)` with
that many top alternatives and the fixed `n_bootstrap`, saving & returning
all resulting figures in a Dict keyed by `m`.
"""
function run_bootstrap_over_alternatives(
        scores_df::DataFrame,
        alt_list::Vector{Int};
        n_bootstrap::Int,
        candidate_cols::Vector{String}=CANDIDATOS_eseb2022,
        demographics::Vector{String}  = ["D02","D10","D12a","Q19_cat"],
        force_include::Vector{String}=String[]
    )

    figs = Dict{Int, Figure}()
    for m in alt_list
        fig = run_bootstrap_analysis(
            scores_df;
            n_bootstrap    = n_bootstrap,
            n_alternatives = m,
            candidate_cols = candidate_cols,
            demographics   = demographics, force_include = force_include
        )
        figs[m] = fig
    end
    return figs
end




"""
    compute_measures_over_alternatives(
        scores_df::DataFrame,
        alt_list::Vector{Int};
        n_bootstrap::Int,
        candidate_cols::Vector{String}=CANDIDATOS_eseb2022,
        demographics::Vector{String}  = ["D02","D10","D12a","Q19_cat"]
    ) → Dict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}

For each `m` in `alt_list`, runs `compute_bootstrap_measures` with that many
top candidates and `n_bootstrap`, returning a Dict mapping `m` → measures-dict.

This avoids any plotting and does only the heavy impute/encode/bootstrap/measure
pipeline. To keep it as performant as possible, we:

1.  Preallocate the output Dict with the right capacity.
2.  Pass through the same `candidate_cols` and `demographics` without mutation.
3.  Avoid any needless intermediate copies inside the loop.
"""
function compute_measures_over_alternatives(
        scores_df::DataFrame,
        alt_list::Vector{Int};
        n_bootstrap::Int,
        candidate_cols::Vector{String}=CANDIDATOS_eseb2022,
        demographics::Vector{String}  = ["D02","D10","D12a","Q19_cat"],
        force_include::Vector{String}=String[]
    )
    # reserve output dict for alt_list size
    results = Dict{Int,Any}(zip(alt_list, repeat([nothing], length(alt_list))))

    for m in alt_list
        # Compute and store measures for top-m candidates
        measures_dict = compute_bootstrap_measures(
            scores_df;
            candidate_cols  = candidate_cols,
            demographics    = demographics,
            n_bootstrap     = n_bootstrap,
            n_alternatives  = m, force_include = force_include
        )
        results[m] = measures_dict
    end

    return results
end




"""
    compute_group_bootstrap_measures(
        scores_df::DataFrame;
        candidate_cols::Vector{String}=CANDIDATOS_eseb2022,
        demographics::Vector{String},
        n_bootstrap::Int,
        n_alternatives::Int
    ) -> Dict{Symbol,Dict{Symbol,Vector{Float64}}}

Runs imputation → encode → bootstrap → group‐based metrics for each
demographic column. Returns a Dict mapping each demographic Symbol to its
stats‐dict (`:C`,`:D`,`:G`).
"""
function compute_group_bootstrap_measures(
        scores_df::DataFrame;
        candidate_cols::Vector{String},
        demographics::Vector{String},
        n_bootstrap::Int,
        n_alternatives::Int, force_include::Vector{String}=String[] 
    )
   
    bt_top = build_bootstrapped_encoded_variants(
        scores_df;
        candidate_cols  = candidate_cols,
        demographics    = demographics,
        n_bootstrap     = n_bootstrap,
        ncandidates     = n_alternatives,
        force_include = force_include
    )

   
    stats = Dict{Symbol,Dict{Symbol,Dict{Symbol, Vector{Float64}}}}()
    pm.@showprogress for demo_str in demographics
        demo_sym = Symbol(demo_str)
        sdict = encoded_bootstrap_group_metrics(bt_top, demo_sym)
        add_G!(sdict)
        stats[demo_sym] = sdict
    end

    return stats
end

function combined_CDG_boxplots(
        stats_vec::Vector;                 # ← NEW
        titles::Vector,
        n_bootstrap::Int,
        n_alternatives::Int,
        variants  = ["zero","random","mice"],
        palette   = Makie.wong_colors(),
        boxwidth  = 0.18,
        ncols::Int = 2,
        figsize   = (1400, 900))

    @assert length(stats_vec) == length(titles) "stats_vec and titles must match"
    K        = length(stats_vec)
    nrows    = ceil(Int, K / ncols)
    measures = [:C, :D, :G]
    mlabels  = ["C (coherence)", "D (divergence)", "G (√C·D)"]
    offsets  = LinRange(-boxwidth, boxwidth, length(measures))
    nv       = length(variants)

    fig = Figure(resolution = figsize)

    # Global title
    Label(fig[0, 1:(ncols)];
          text     = "B = $n_bootstrap    •    m = $n_alternatives",
          fontsize = 18,
          halign   = :center)

    # Axes grid
    axes = Dict{Int,Axis}()
    for p in 1:K
        row, col = fldmod1(p, ncols)     # 1-based (row, col)
        axes[p] = Axis(fig[row, col];
                       title  = titles[p],
                       xlabel = "imputation variant",
                       ylabel = "value")
    end

    legend_handles = BoxPlot[]; legend_labels = String[]

    for (p, stats) in pairs(stats_vec)                # panel loop
        ax = axes[p]
        for (j, meas) in enumerate(measures)          # C / D / G
            xs = Float32[]; ys = Float32[]
            for (i, v) in enumerate(variants)         # zero / random / mice
                vals = stats[Symbol(v)][meas]
                append!(xs, fill(i + offsets[j], length(vals)))
                append!(ys, Float32.(vals))
            end
            col = palette[(j-1) % length(palette) + 1]
            bp  = boxplot!(ax, xs, ys; color = col, width = boxwidth*0.9)

            if p == 1                                 # capture legend once
                push!(legend_handles, bp)
                push!(legend_labels, mlabels[j])
            end
        end
        ax.xticks = (1:nv, variants)
    end

    # Legend in an extra column on the right
    Legend(fig[1:nrows, ncols+1], legend_handles, legend_labels)

    return fig
end


function run_group_bootstrap_analysis(
        scores_df::DataFrame;
        n_bootstrap::Int,
        n_alternatives::Int,
        candidate_cols::Vector{String}=CANDIDATOS_eseb2022,
        demographics::Vector{String},
        force_include::Vector{String}=String[]
    )

    # 1 ─ stats for every demographic in *this* order
    stats = compute_group_bootstrap_measures(
        scores_df;
        candidate_cols  = candidate_cols,
        demographics    = demographics,
        n_bootstrap     = n_bootstrap,
        n_alternatives  = n_alternatives,
        force_include = force_include)

    # 2 ─ collect to vectors expected by the new combined_CDG_boxplots
    stats_vec = [stats[Symbol(d)] for d in demographics]
    titles    = demographics                      # or prettified labels

    fig = combined_CDG_boxplots(
              stats_vec;                          # ← new positional arg
              titles         = titles,
              n_bootstrap    = n_bootstrap,
              n_alternatives = n_alternatives)

    fname = "exploratory/imgs/bootstrapTop$(n_alternatives)_CDG_B$(n_bootstrap).png"
    save(fname, fig; px_per_unit = 2)
    println("Saved group-based plot → ", fname)

    return fig 
end


"""
    run_group_over_alternatives(
        scores_df::DataFrame,
        alt_list::Vector{Int};
        n_bootstrap::Int,
        candidate_cols::Vector{String}=CANDIDATOS_eseb2022,
        demographics::Vector{String}
    ) -> Dict{Int, Figure}

Loops `run_group_bootstrap_analysis` over each `m` in `alt_list`, returning
a Dict mapping `m` → the saved Figure.
"""
function run_group_over_alternatives(
        scores_df::DataFrame,
        alt_list::Vector{Int};
        n_bootstrap::Int,
        candidate_cols::Vector{String}=CANDIDATOS_eseb2022,
        demographics::Vector{String},
        force_include::Vector{String}=String[] 
    )

    figs = Dict{Int,Figure}()
    pm.@showprogress for m in alt_list
        figs[m] = run_group_bootstrap_analysis(
            scores_df;
            n_bootstrap     = n_bootstrap,
            n_alternatives  = m,
            candidate_cols  = candidate_cols,
            demographics    = demographics,
            force_include = force_include
        )
    end
    return figs
end


function plot_group_over_alternatives(
        scores_df,
        demo,
        alt_list::Vector{Int};
        candidate_cols::Vector{String}=CANDIDATOS_eseb2022,
        n_bootstrap::Int              = 500,
        boxwidth::Real                = 0.18,
        palette                       = Makie.wong_colors(),
        figsize                       = (850, 550),
        force_include::Vector{String}=String[] )

    figs = Dict{Int,Figure}()

    variants       = ["zero","random","mice"]
    measures       = [:C, :D, :G]
    measure_labels = ["C", "D", "G"]
    n_variants     = length(variants)

   pm.@showprogress for m in alt_list
        # 1 ─ bootstrap → group stats → add G
        bt    = build_bootstrapped_encoded_variants(scores_df;
                   candidate_cols = candidate_cols,
                   demographics   = [String(demo)],
                   n_bootstrap    = n_bootstrap,
                   ncandidates    = m, force_include = force_include)

        stats = encoded_bootstrap_group_metrics(bt, demo)
        add_G!(stats)                     # adds :G vectors in-place

        # 2 ─ figure: 1 row × 2 columns  [ axis | legend ]
        fig  = Figure(resolution = figsize)
        grid = fig[1, 1] = GridLayout(cols = (Relative(0.75), Relative(0.25)))

        # title over both columns
        Label(grid[0, 1:2],
              "$(demo) — B=$(n_bootstrap) • m=$(m)",
              fontsize = 16, halign = :left)

        # main axis
        ax = Axis(grid[1, 1],
                  xlabel = "imputation variant",
                  ylabel = "value")

        offsets = LinRange(-boxwidth, boxwidth, length(measures))

        legend_handles = Any[]; legend_labels = String[]

        for (j, meas) in enumerate(measures)
            xs = Float32[]; ys = Float32[]
            for (i, var) in enumerate(variants)
                vals = stats[Symbol(var)][meas]
                append!(xs, fill(i + offsets[j], length(vals)))
                append!(ys, Float32.(vals))
            end
            col = palette[(j-1) % length(palette) + 1]
            h   = boxplot!(ax, xs, ys; color = col, width = boxwidth*0.9)
            push!(legend_handles, h)
            push!(legend_labels, measure_labels[j])
        end

        ax.xticks = (1:n_variants, variants)

        # legend in the right-hand column
        Legend(grid[1, 2], legend_handles, legend_labels)

        figs[m] = fig
    end

    return figs
end



function plot_margin_graph_from_bootstrap(bt_variant::Dict; title="", digits=2)
    stats = summarize_margins(margins_over_bootstrap(bt_variant))
    return _margin_graph_py(stats; title, digits)
end

function plot_margin_graph_from_table(df::DataFrame; title="", digits=2)
    # turn the DataFrame into Dict{Tuple,Tuple} like the other path
    stats = Dict{Tuple{Symbol,Symbol},Tuple{Float64,Float64}}(
        (Symbol(r.cand1), Symbol(r.cand2)) => (r.mean_margin, r.sd_margin)
        for r in eachrow(df))
    return _margin_graph_py(stats; title, digits)
end



# ── Dict → square DataFrame ───────────────────────────────────────────────
function margin_matrix(stats::Dict{Tuple{Symbol,Symbol},Tuple{Float64,Float64}};
                       digits::Int = 3)

    cands = sort!(unique(vcat(first.(keys(stats)), last.(keys(stats)))))
    K     = length(cands)
    idx   = Dict(c => i for (i,c) in enumerate(cands))

    mat = fill(0.0, K, K)
    for ((a,b),(μ,_)) in stats
        i, j = idx[a], idx[b]
        μ    = round(μ; digits=digits)
        mat[i,j] =  μ
        mat[j,i] = -μ
    end

    df = DataFrame(mat, Symbol.(cands))
    insertcols!(df, 1, :cand => cands)
    return df
end


# ── long DataFrame → square DataFrame ─────────────────────────────────────
function margin_matrix(df_long::DataFrame; digits::Int = 3)
    cands = sort!(unique(vcat(df_long.cand1, df_long.cand2)))
    K     = length(cands)
    idx   = Dict(c => i for (i,c) in enumerate(cands))

    mat = fill(0.0, K, K)
    for r in eachrow(df_long)
        i, j = idx[r.cand1], idx[r.cand2]
        μ    = round(r.mean_margin; digits=digits)
        mat[i,j] =  100*μ
        mat[j,i] = 100*-μ
    end

    df = DataFrame(mat, Symbol.(cands))
    insertcols!(df, 1, :cand => cands)
    return df
end



function compute_stats_over_m(
    scores_df;
    candidate_cols,
    demographics,
    m_values    = 2:7,
    n_bootstrap = 750,
    force_include
)

    stats_m = Dict{Int,Dict{Symbol,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}}()

    for m in m_values
        stats = compute_group_bootstrap_measures(
            scores_df;
            candidate_cols  = candidate_cols,
            demographics    = demographics,
            n_bootstrap     = n_bootstrap,
            n_alternatives  = m,
            force_include = force_include
        )
        add_G!.(values(stats))             # mutate in place → adds :G
        stats_m[m] = stats
    end
    return stats_m
end
