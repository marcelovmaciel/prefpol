# ─────────────────────────────────────────────────────────────────────────────
# Scenario inside a TOML file
# ─────────────────────────────────────────────────────────────────────────────
struct Scenario
    name       :: String
    candidates :: Vector{String}
end

# ─────────────────────────────────────────────────────────────────────────────
# FULL election specification (= everything in the TOML)
# ─────────────────────────────────────────────────────────────────────────────
struct ElectionConfig
    year            :: Int
    data_loader     :: String
    data_file       :: String   
    max_candidates  :: Int
    m_values_range  :: Vector{Int}

    n_bootstrap     :: Int          # default; can be overridden later
    n_alternatives  :: Int
    rng_seed        :: Int

    candidates      :: Vector{String}
    demographics    :: Vector{String}
    scenarios       :: Vector{Scenario}    # list of Scenario structs
end

# ─────────────────────────────────────────────────────────────────────────────
# RUN-TIME configuration for *one* bootstrap run
# ─────────────────────────────────────────────────────────────────────────────
struct BootstrapConfig
    election   :: ElectionConfig   # full context
    scenario   :: Scenario         # which scenario
    n_bootstrap:: Int              # number of replicates for *this* run
    rng_seed   :: Int              # seed for *this* run
end

function load_election_cfg(path::AbstractString)::ElectionConfig
    t = TOML.parsefile(path)
    proj = dirname(Pkg.project().path)       # project root
    rawfile = isabspath(t["data_file"]) ? t["data_file"] :
              joinpath(proj, t["data_file"])

    scen_vec = [Scenario(s["name"], Vector{String}(s["candidates"]))
                for s in t["forced_scenarios"]]

    ElectionConfig(
        t["year"], t["data_loader"], rawfile,
        t["max_candidates"], Vector{Int}(t["m_values_range"]),
        t["n_bootstrap"], t["n_alternatives"],
        t["rng_seed"],
        Vector{String}(t["candidates"]),
        Vector{String}(t["demographics"]),
        scen_vec,
    )
end



import Base: show



const IND = "    "                       # 4-space indent

# ────────────────────────── helpers ──────────────────────────
_pp(io, key, val, lvl) = println(io, repeat(IND, lvl), key, " = ", val)
function _pp_vec(io, key, vec, lvl; max=8)
    head = first(vec, max)
    tail = length(vec) > max ? " … ("*string(length(vec))*")" : ""
    println(io, repeat(IND, lvl), key, " = ", head, tail)
end

# ────────────────────────── Scenario ─────────────────────────
function show(io::IO, ::MIME"text/plain", s::Scenario; kwargs...)
    print(io, "Scenario(\"", s.name, "\", ", s.candidates, ")")
end

# ─────────────────────── ElectionConfig ──────────────────────
function show(io::IO, ::MIME"text/plain", ec::ElectionConfig; kwargs...)
    println(io, "ElectionConfig(")
    _pp(io, "year",             ec.year, 1)
    _pp(io, "data_loader",      ec.data_loader, 1)
    _pp(io, "data_file",        ec.data_file, 1)
    _pp(io, "max_candidates",   ec.max_candidates, 1)
    _pp(io, "m_values_range",   ec.m_values_range, 1)
    _pp(io, "n_bootstrap(def)", ec.n_bootstrap, 1)
    _pp(io, "n_alternatives",   ec.n_alternatives, 1)
    _pp(io, "rng_seed(def)",    ec.rng_seed, 1)
    _pp_vec(io, "candidates",   ec.candidates, 1)
    _pp(io, "demographics",     ec.demographics, 1)
    println(io, IND, "scenarios = [")
    for sc in ec.scenarios
        println(io, repeat(IND,2), sc)           # uses Scenario show
    end
    println(io, IND, "]")
    print(io, ")")
end

# ─────────────────────── BootstrapConfig ─────────────────────
function show(io::IO, ::MIME"text/plain", bc::BootstrapConfig; kwargs...)
    println(io, "BootstrapConfig(")
    _pp(io, "year",        bc.election.year, 1)
    _pp(io, "scenario",    bc.scenario.name, 1)
    _pp(io, "forced_set",  bc.scenario.candidates, 1)
    _pp(io, "n_bootstrap", bc.n_bootstrap, 1)
    _pp(io, "rng_seed",    bc.rng_seed, 1)
    _pp(io, "m_values",    bc.election.m_values_range, 1)
    _pp(io, "data_file",   bc.election.data_file, 1)
    print(io, IND, "election = ")
    show(io, MIME"text/plain"(), bc.election; kwargs...)  # pass through
    println(io)
    print(io, ")")
end

"""
    make_bootstrap_cfg(elec, scenario_name; n_bootstrap=elec.n_bootstrap,
                       rng_seed = elec.rng_seed)

Creates the concrete run-time configuration for **one** bootstrap/analysis job.
"""
function make_bootstrap_cfg(elec::ElectionConfig, scenario_name::AbstractString;
                            n_bootstrap::Int = elec.n_bootstrap,
                            rng_seed   ::Int = elec.rng_seed)

    scen = findfirst(s -> s.name == scenario_name, elec.scenarios)
    scen === nothing && error("Scenario $scenario_name not found for year $(elec.year)")

    BootstrapConfig(elec, elec.scenarios[scen], n_bootstrap, rng_seed)
end




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


function load_election_data(cfg::ElectionConfig)
    loader_sym = Symbol(cfg.data_loader)

    # resolve in the current module’s namespace
    if !isdefined(@__MODULE__, loader_sym)
        throw(ArgumentError("data_loader ‘$(cfg.data_loader)’ not found in module $(nameof(@__MODULE__))"))
    end

    loader_fun = getfield(@__MODULE__, loader_sym)
    return loader_fun(cfg.data_file; candidates = cfg.candidates)
end

# ─────────────────────────────────────────────────────────────────────────────
# Candidate set for ONE BootstrapConfig, using your existing helper
# ─────────────────────────────────────────────────────────────────────────────
"""
    candidate_set_for(bc::BootstrapConfig, df::DataFrame) -> Vector{String}

Uses `compute_candidate_set` with
* all survey candidates (`bc.election.candidates`);
* `m = bc.election.n_alternatives`;
* `force_include = bc.scenario.candidates`.
"""
function candidate_set_for(bc::BootstrapConfig, df::DataFrame)
    compute_candidate_set(df;
        candidate_cols = bc.election.candidates,
        m              = bc.election.n_alternatives,
        force_include  = bc.scenario.candidates)
end


function build_bootstrapped_encoded_variants(bc::BootstrapConfig)
df = load_election_data(bc.election)
cset = candidate_set_for(bc, df)   
return build_bootstrapped_encoded_variants(
           df;
           candidate_cols = bc.election.candidates,      # or bc.election.candidates
           demographics   = bc.election.demographics,
           n_bootstrap    = bc.n_bootstrap,
           ncandidates    = length(cset),
           force_include  = cset,
       )
end

# ════════════════════════════════════════════════════════════════════════════
# 0 · prerequisites
# ════════════════════════════════════════════════════════════════════════════

const INT_DIR = "intermediate_data"
mkpath(INT_DIR)

# ════════════════════════════════════════════════════════════════════════════
# 1 · lightweight container  +  deterministic filename helpers
# ════════════════════════════════════════════════════════════════════════════
struct BootstrapReplicates
    cfg       :: BootstrapConfig
    variants  :: Dict{Symbol,Any}          # :mice → Dict with :dfs, …
    timestamp :: DateTime
end

bootfile(bc::BootstrapConfig) =
    joinpath(INT_DIR, "boot_$(bc.election.year)_$(bc.scenario.name).jld2")

function save_bootstraps(br::BootstrapReplicates)
    @save bootfile(br.cfg) br
    br
end

function load_bootstraps(bc::BootstrapConfig)::BootstrapReplicates
    @load bootfile(bc) br
    return br
end



# ─────────────────────────────────────────────────────────────────────────────
#  Helper: discover *.toml files and map to their YEAR key
# ─────────────────────────────────────────────────────────────────────────────
function _discover_tomls(cfgdir::String)
    pairs = Dict{Int,String}()
    for f in filter(x -> endswith(x, ".toml"), readdir(cfgdir, join=true))
        year = try
            # safer: read TOML, fetch its year key
            TOML.parsefile(f)["year"]
        catch
            @warn "Skipping $f (no readable ‘year’ key)"; continue
        end
        pairs[year] = f
    end
    return pairs
end

# ─────────────────────────────────────────────────────────────────────────────
#  Batch driver
# ─────────────────────────────────────────────────────────────────────────────
"""
    build_boootstraped_scenarios(; years   = nothing,
                                 which   = nothing,   # Dict(year => [scenario …])
                                 n_boot  = nothing,   # Int | nothing
                                 cache   = true,
                                 cfgdir  = "config")

Run all requested `(year, scenario)` combinations, caching bootstraps
in `intermediate_data/`.

Returns `Dict{Int,Dict{String,BootstrapReplicates}}`.
"""
function build_boootstraped_scenarios(; years = nothing,
                                      which = nothing,
                                      n_boot = nothing,
                                      cache = true,
                                      cfgdir = "config")

    toml_by_year = _discover_tomls(cfgdir)

    years_vec = years === nothing      ? sort(collect(keys(toml_by_year))) :
                 isa(years, Integer)   ? [years] :
                 years

    out = Dict{Int,Dict{String,BootstrapReplicates}}()

    for yr in years_vec
        haskey(toml_by_year, yr) || (@warn "No config for year $yr"; continue)

        elec = load_election_cfg(toml_by_year[yr])
        raw  = load_election_data(elec)

        # Determine which scenarios to run
        scen_names = which === nothing           ? [s.name for s in elec.scenarios] :
                     haskey(which, yr)           ? which[yr] :
                     ( @warn "No scenario list for year $yr in `which=`"; String[] )

        res = Dict{String,BootstrapReplicates}()

        for sn in scen_names
            bc = make_bootstrap_cfg(
                    elec, sn;
                    n_bootstrap = n_boot === nothing ? elec.n_bootstrap : n_boot)

            # choose cached or fresh
            br = cache && isfile(bootfile(bc)) ?
                     ( @info "• $yr/$sn  loading cache"; load_bootstraps(bc) ) :
                     begin
                         @info "• $yr/$sn  generating"
                         reps = build_bootstrapped_encoded_variants(bc)              # ← your 1-arg method
                         save_bootstraps(BootstrapReplicates(bc, reps, now()))
                     end

            res[sn] = br
        end
        out[yr] = res
    end
    return out
end


# ================ Below here is what I had written in pipeline ================
# ════════════════════════════════════════════════════════════════════════════


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




# ================ Below here is what claude had written ================
# ════════════════════════════════════════════════════════════════════════════












struct AnalysisResults
    config::BootstrapConfig
    global_measures::Any
    group_measures::Any
    bootstraps::Dict{Symbol,Vector{DataFrame}}      # NEW field

    function AnalysisResults(config::BootstrapConfig,
                             global_measures,
                             group_measures,
                             bootstraps::Dict{Symbol,Vector{DataFrame}})
        new(config, global_measures, group_measures, bootstraps)
    end
end
# Plot Results Structure
struct PlotResults
    analysis_results::AnalysisResults
    global_plot::Any
    group_plot::Any
    
    function PlotResults(analysis_results::AnalysisResults, global_plot, group_plot)
        new(analysis_results, global_plot, group_plot)
    end
end

# Data Loading Module
struct DataLoader
    config::Dict
    
    function DataLoader()
        new(ESEBConfig.ELECTIONS)
    end
end


function compute_candidate_sets(loader::DataLoader, data, year::Int)
    election_config = loader.config[year]
    candidate_sets = Dict{String, Vector{String}}()
    
    for scenario in election_config[:forced_scenarios]
        scenario_name = scenario[:name]
        forced_candidates = scenario[:candidates]
        
        if isempty(forced_candidates)
            candidate_set = pp.compute_candidate_set(
                data;
                candidate_cols = election_config[:candidates],
                m = election_config[:max_candidates]
            )
        else
            candidate_set = pp.compute_candidate_set(
                data;
                candidate_cols = election_config[:candidates],
                m = election_config[:max_candidates],
                force_include = forced_candidates
            )
        end
        
        candidate_sets[scenario_name] = candidate_set
    end
    
    return candidate_sets
end

# Analysis Engine - Only computes measures
struct MeasureCalculator
    data_loader::DataLoader
    
    function MeasureCalculator()
        new(DataLoader())
    end
end

function compute_global_measures(calculator::MeasureCalculator, data, config::BootstrapConfig)
    @info "Computing global measures for $(config.year) $(config.scenario_name)..."
    
    kwargs = Dict(
        :n_bootstrap => config.n_bootstrap,
        :candidate_cols => config.candidate_cols,
        :demographics => config.demographics
    )
    
    if !isempty(config.forced_candidates)
        kwargs[:force_include] = config.forced_candidates
    end
    
    return pp.compute_measures_over_alternatives(data, config.m_values; kwargs...)
end

function compute_group_measures(calculator::MeasureCalculator, data, config::BootstrapConfig)
    @info "Computing group measures for $(config.year) $(config.scenario_name)..."
    
    kwargs = Dict(
        :candidate_cols => config.candidate_cols,
        :demographics => config.demographics,
        :m_values => config.m_values,
        :n_bootstrap => config.n_bootstrap,
        :force_include => config.forced_candidates
    )
    
    return pp.compute_stats_over_m(data; kwargs...)
end



"""
    generate_bootstraps(df, cfg) → Dict{Symbol,Vector{DataFrame}}

• Applies the three imputation variants (zero, random, mice)  
• Draws `cfg.n_bootstrap` weighted replicates for each variant  
• Returns the `Dict(:zero=>[..], :random=>[..], :mice=>[..])`
"""
function generate_bootstraps(df::DataFrame, cfg::BootstrapConfig)
    imps = pp.imputation_variants(df,               # ← the imputation you already wrote
                               cfg.candidate_cols,
                               cfg.demographics;
                               most_known_candidates = cfg.candidate_set)

    # survey weight column is usually called :peso or :weight – adapt if needed
    w = hasproperty(df, :peso) ? Vector{Float64}(df.peso) :
        hasproperty(df, :weight) ? Vector{Float64}(df.weight) :
        ones(Float64, nrow(df))            # fallback: equal weights

    return pp.bootstrap_variants(imps; weights = w, B = cfg.n_bootstrap)
end


function analyze_scenario(calculator::MeasureCalculator,
                          data::DataFrame,
                          cfg::BootstrapConfig)

    @info "Analyzing scenario: $(cfg.year) – $(cfg.scenario_name)"

    # ─── ❶ generate & store replicates ─────────────────────────
    boots = generate_bootstraps(data, cfg)

    # ─── ❷ your existing point-estimates ───────────────────────
    global_meas = compute_global_measures(calculator, data, cfg)
    group_meas  = compute_group_measures(calculator,  data, cfg)

    # ─── ❸ bundle everything ──────────────────────────────────
    return AnalysisResults(cfg, global_meas, group_meas, boots)
end



# Plotting Engine - Separate from calculation
struct PlotGenerator
    output_dir::String
    
    function PlotGenerator(output_dir::String = "exploratory/imgs")
        new(output_dir)
    end
end

function create_global_plot(generator::PlotGenerator, results::AnalysisResults; variants = ["zero", "random", "mice"])
    config = results.config
    candidate_label = pp.describe_candidate_set(config.candidate_set)
    
    return pp.lines_alt_by_variant(
        results.global_measures, 
        candidate_label = candidate_label, 
        year = config.year,
        variants = variants
    )
end

function create_group_plot(generator::PlotGenerator, results::AnalysisResults; variants = ["zero", "random", "mice"])
    config = results.config
    candidate_label = pp.describe_candidate_set(config.candidate_set)
    maxcols = config.year == 2006 ? 1 : 2
    
    return pp.lines_group_measures_over_m(
        results.group_measures,
        demographics = Symbol.(config.demographics),
        candidate_label = candidate_label,
        year = config.year,
        maxcols = maxcols,
        variants = variants
    )
end

function create_plots(generator::PlotGenerator, results::AnalysisResults; variants = ["zero", "random", "mice"])
    @info "Creating plots for $(results.config.year) $(results.config.scenario_name)..."
    
    global_plot = create_global_plot(generator, results, variants = variants)
    group_plot = create_group_plot(generator, results, variants = variants)
    
    return PlotResults(results, global_plot, group_plot)
end

function save_plots(generator::PlotGenerator, plot_results::PlotResults)
    config = plot_results.analysis_results.config
    bootstrap_suffix = "$(config.n_bootstrap)bts"
    
    global_filename = "alt_by_variant_lines$(config.year)_$(bootstrap_suffix)_$(config.scenario_name).png"
    group_filename = "cgd_lines$(config.year)_$(bootstrap_suffix)_$(config.scenario_name).png"
    
    global_path = joinpath(generator.output_dir, global_filename)
    group_path = joinpath(generator.output_dir, group_filename)
    
    pp.save(global_path, plot_results.global_plot; px_per_unit=3)
    pp.save(group_path, plot_results.group_plot; px_per_unit=3)
    
    @info "Saved plots: $global_filename, $group_filename"
end

# Configuration Builder Functions
function create_bootstrap_configs(year::Int, n_bootstrap::Int = 1500)
    loader = DataLoader()
    data = load_election_data(loader, year)
    candidate_sets = compute_candidate_sets(loader, data, year)
    
    configs = BootstrapConfig[]
    election_config = ESEBConfig.ELECTIONS[year]
    
    for scenario in election_config[:forced_scenarios]
        scenario_name = scenario[:name]
        candidate_set = candidate_sets[scenario_name]
        forced_candidates = scenario[:candidates]
        
        config = BootstrapConfig(
            year = year,
            scenario_name = scenario_name,
            candidate_set = candidate_set,
            forced_candidates = forced_candidates,
            n_bootstrap = n_bootstrap,
            m_values = collect(election_config[:m_values_range]),
            demographics = election_config[:demographics],
            candidate_cols = election_config[:candidates],
            rng_seed = election_config[:rng_seed])
        
        push!(configs, config)
    end
    
    return configs, data
end

# High-level Workflow Functions
function run_analysis_pipeline(configs::Vector{BootstrapConfig}, data; 
                             create_plots_flag::Bool = true, 
                             save_plots_flag::Bool = true,
                             plot_variants = ["mice"])
    calculator = MeasureCalculator()
    plot_generator = PlotGenerator()
    
    analysis_results = AnalysisResults[]
    plot_results = PlotResults[]
    
    # Calculate measures for all configurations
    @info "Computing measures for $(length(configs)) configurations..."
    for config in configs
        result = analyze_scenario(calculator, data, config)
        push!(analysis_results, result)
    end
    
    # Create plots if requested
    if create_plots_flag
        @info "Creating plots..."
        for result in analysis_results
            plot_result = create_plots(plot_generator, result, variants = plot_variants)
            push!(plot_results, plot_result)
            
            if save_plots_flag
                save_plots(plot_generator, plot_result)
            end
        end
    end
    
    return analysis_results, plot_results
end

function analyze_year(year::Int; n_bootstrap::Int = 1500, create_plots_flag::Bool = true, save_plots_flag::Bool = true)
    @info "Starting analysis for year $year with $n_bootstrap bootstrap samples"
    
    configs, data = create_bootstrap_configs(year, n_bootstrap)
    return run_analysis_pipeline(configs, data; 
                               create_plots_flag = create_plots_flag, 
                               save_plots_flag = save_plots_flag)
end

function analyze_multiple_years(years::Vector{Int}; n_bootstrap::Int = 1500, create_plots_flag::Bool = true, save_plots_flag::Bool = true)
    all_analysis_results = Dict{Int, Vector{AnalysisResults}}()
    all_plot_results = Dict{Int, Vector{PlotResults}}()
    
    for year in years
        analysis_results, plot_results = analyze_year(year; 
                                                    n_bootstrap = n_bootstrap,
                                                    create_plots_flag = create_plots_flag,
                                                    save_plots_flag = save_plots_flag)
        all_analysis_results[year] = analysis_results
        all_plot_results[year] = plot_results
    end
    
    return all_analysis_results, all_plot_results
end

# Convenience functions for different use cases
function run_development_analysis(years::Vector{Int} = [2022])
    return analyze_multiple_years(years; n_bootstrap = 10, create_plots_flag = true, save_plots_flag = false)
end

function run_production_analysis(years::Vector{Int} = [2006, 2018, 2022])
    return analyze_multiple_years(years; n_bootstrap = 1500, create_plots_flag = true, save_plots_flag = true)
end

function run_quick_test(year::Int = 2022)
    return analyze_year(year; n_bootstrap = 5, create_plots_flag = false, save_plots_flag = false)
end

# Analysis-only functions (no plotting)
function compute_measures_only(years::Vector{Int}; n_bootstrap::Int = 1500)
    return analyze_multiple_years(years; n_bootstrap = n_bootstrap, create_plots_flag = false, save_plots_flag = false)
end

# Plotting-only functions (from existing results)
function create_plots_from_results(analysis_results::Vector{AnalysisResults}; 
                                 variants = ["mice"], 
                                 save_plots_flag::Bool = true)
    plot_generator = PlotGenerator()
    plot_results = PlotResults[]
    
    for result in analysis_results
        plot_result = create_plots(plot_generator, result, variants = variants)
        push!(plot_results, plot_result)
        
        if save_plots_flag
            save_plots(plot_generator, plot_result)
        end
    end
    
    return plot_results
end

# Main execution function
function main()
    # Example usage:
    
    # 1. For development - fast execution, no saving
    # analysis_results, plot_results = run_development_analysis([2022])
    
    # 2. For production - full analysis with plots
    # analysis_results, plot_results = run_production_analysis()
    
    # 3. Compute measures only (for later plotting)
    # analysis_results, _ = compute_measures_only([2022])
    
    # 4. Create plots from existing analysis results
    # plot_results = create_plots_from_results(analysis_results[2022])
    
    # Default: production analysis
    return run_production_analysis()
end




function compute_bootstraps_only(year::Int; n_bootstrap::Int = 1500)
    # ❶ Data + configs
    configs, data = create_bootstrap_configs(year, n_bootstrap)
    results = Dict{String,Dict}()

    # ❷ For each scenario: generate replicates, optionally save them
    for cfg in configs
        br = generate_bootstraps(data, cfg)                    # ← fast
        results[cfg.scenario_name] = br
        # uncomment next line if you also want an on-disk copy
        # save_bootstraps("exploratory/imgs/bootstraps_$(year)_$(cfg.scenario_name).jld2", br)
    end
    return results
end



