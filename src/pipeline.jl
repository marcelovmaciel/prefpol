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




function load_election_data(cfg::ElectionConfig)
    loader_sym = Symbol(cfg.data_loader)

    # resolve in the current module’s namespace
    if !isdefined(@__MODULE__, loader_sym)
        throw(ArgumentError("data_loader ‘$(cfg.data_loader)’ not found in module $(nameof(@__MODULE__))"))
    end

    loader_fun = getfield(@__MODULE__, loader_sym)
    return loader_fun(cfg.data_file; candidates = cfg.candidates)
end




function weighted_bootstrap(ecfg::ElectionConfig)
    df = load_election_data(ecfg)
    weights = df.peso 
    B = ecfg.n_bootstrap
    # slice df to only candidates and demographics from ecfg
    candidates = ecfg.candidates
    demographics = ecfg.demographics
    df = select(df, candidates..., demographics...)
    bts = weighted_bootstrap(df, weights, B)
    return bts
end




const INT_DIR = "intermediate_data"
mkpath(INT_DIR)

# ════════════════════════════════════════════════════════════════════════════
# 1 · lightweight container  +  deterministic filename helpers
# ════════════════════════════════════════════════════════════════════════════
struct BootstrapReplicates
    cfg       :: ElectionConfig
    timestamp :: DateTime
end

bootfile(bc::ElectionConfig) =
    joinpath(INT_DIR, "boot_$(bc..year).jld2")

function save_bootstraps(br::BootstrapReplicates)
    @save bootfile(br.cfg) br
    br
end

function load_bootstraps(bc::ElectionConfig)::BootstrapReplicates
    @load bootfile(bc) br
    return br
end



"""
    save_bootstrap(cfg; dir = INT_DIR, overwrite = false, quiet = false)
        → NamedTuple{(:path,:data,:cached), ...}

Ensure a weighted-bootstrap exists for `cfg`:

  • If `dir/boot_YEAR.jld2` is missing *or* `overwrite=true`, build the
    bootstrap with `weighted_bootstrap(cfg)` and write it to disk.

  • Otherwise reuse the cached file, **loading the replicates** so that
    `.data` is never `nothing`.

Returned fields
---------------
| field   | meaning                                   |
|---------|-------------------------------------------|
| `path`  | full path to the `.jld2` file             |
| `data`  | the `reps` object (always in memory)      |
| `cached`| `true` if we reused an existing file      |
"""
function save_bootstrap(cfg::ElectionConfig;
                        dir::AbstractString = INT_DIR,
                        overwrite::Bool = false,
                        quiet::Bool = false)

    path = joinpath(dir, "boot_$(cfg.year).jld2")

    # ------------------ cache hit ------------------
    if !overwrite && isfile(path)
        !quiet && @warn "Reusing cached bootstrap at $(path); loading into memory"
        reps = nothing
        @load path reps                           # brings `reps` back
        return (path = path, data = reps, cached = true)
    end

    # ------------------ (re)build ------------------
    reps = weighted_bootstrap(cfg)                # heavy call
    @save path reps cfg
    !quiet && @info "Saved bootstrap for year $(cfg.year) → $(path)"
    return (path = path, data = reps, cached = false)
end


# ————————————————————————————————————————————————————————————————
#  2.  Batch driver: all (or selected) years in one call
# ————————————————————————————————————————————————————————————————
"""
    save_all_bootstraps(; years = nothing,
                            cfgdir = "config",
                            overwrite = false) -> Dict{Int,String}

Iterate over every `*.toml` in `cfgdir`; for each year that matches
`years` (or *all* years if `years === nothing`) build & save the bootstrap.

Returns a dictionary `year ⇒ saved_filepath`.
"""
function save_all_bootstraps(; years = nothing,
                             cfgdir::AbstractString = "config",
                             overwrite::Bool = false)
    # discover configs on disk
    toml_files = filter(p -> endswith(p, ".toml"), readdir(cfgdir; join=true))
    isempty(toml_files) && error("No TOML files found in $(cfgdir)")

    wanted = years === nothing        ? nothing           :
             isa(years, Integer)      ? Set([years])      :
             Set(years)

    saved = Dict{Int,String}()

    for f in sort(toml_files)
        cfg = load_election_cfg(f)
        (wanted !== nothing && !(cfg.year in wanted)) && continue
        @info "Processing year $(cfg.year) from $(f)"
        saved[cfg.year] = save_bootstrap(cfg; overwrite).path
    end
    return saved
end






"""
    load_all_bootstraps(; years = nothing,
                           dir   = INT_DIR,
                           quiet = false)
        → OrderedDict{Int,NamedTuple}

Read every `boot_YYYY.jld2` in `dir` (or just the chosen `years`)
and return them in a year-sorted `OrderedDict`.

Each value is a `NamedTuple` with

| field   | meaning                     |
|---------|-----------------------------|
| `data`  | the bootstrap replicates    |
| `cfg`   | the `ElectionConfig` object |
| `path`  | full path to the file       |
"""
function load_all_bootstraps(; years   = nothing,
                             dir::AbstractString = INT_DIR,
                             quiet::Bool = false)

    paths = filter(p -> occursin(r"boot_\d+\.jld2$", p),
                   readdir(dir; join = true))

    isempty(paths) && error("No bootstrap files found in $(dir)")

    selected = years === nothing       ? nothing :
               isa(years,Integer)      ? Set([years]) :
               Set(years)

    out = OrderedCollections.OrderedDict{Int,NamedTuple}()

    for f in sort(paths)                       # alphabetical = chronological
        yr = parse(Int, match(r"boot_(\d{4})\.jld2", basename(f)).captures[1])
        (selected !== nothing && !(yr in selected)) && continue

        reps = cfg = nothing
        @load f reps cfg

        !quiet && @info "Loaded bootstrap $(yr)  ←  $(f)"

        out[yr] = (data = reps, cfg = cfg, path = f)
    end
    return out
end


function add_imputation_variants_to_bts(bt::NamedTuple;
most_known_candidates::Vector{String}=String[])

reps          = bt.data                       # Vector{DataFrame}
cfg           = bt.cfg
B             = length(reps)

variants = Dict{Symbol, Vector{DataFrame}}(
    :zero   => Vector{DataFrame}(undef, B),
    :random => Vector{DataFrame}(undef, B),
    :mice   => Vector{DataFrame}(undef, B),
)

for (i, df) in enumerate(reps)
    imp = imputation_variants(df,
                              cfg.candidates,
                              cfg.demographics;
                              most_known_candidates)
    variants[:zero][i]   = imp.zero
    variants[:random][i] = imp.random
    variants[:mice][i]   = imp.mice
end

return (data = variants, cfg = cfg, path = bt.path)
end



function impute_and_save(bt::NamedTuple;
                         dir::AbstractString = INT_DIR,
                         overwrite::Bool     = false,
                         most_known_candidates::Vector{String} = String[])

    year  = bt.cfg.year
    path  = joinpath(dir, "boot_imp_$(year).jld2")

    if !overwrite && isfile(path)
        @info "Using cached imputed bootstrap at $(path)"
        return path
    end

    imp   = add_imputation_variants_to_bts(bt; most_known_candidates)
    @save path imp
    @info "Saved imputed bootstrap for year $(year) → $(path)"

    # explicit cleanup
    imp = nothing
    GC.gc()

    return path
end


#= # ————————————————————————————————————————————————————————————————
#  2.  Batch driver: iterate over *all* (or selected) years
# ————————————————————————————————————————————————————————————————
"""
    impute_all_bootstraps(; years         = nothing,
                             base_dir      = INT_DIR,
                             imp_dir       = INT_DIR,
                             overwrite     = false,
                             most_known_candidates = String[])
        → OrderedDict{Int,String}

Load every bootstrap in `base_dir`, run imputations, write each year’s
result to `imp_dir`, and immediately release memory.

`years` may be:
  * `nothing`  – all years present on disk;
  * an `Int`   – a single year;
  * `Vector{Int}` – a specific set of years.

Returns an `OrderedDict year ⇒ saved_path`, sorted chronologically.
"""
function impute_all_bootstraps(; years    = nothing,
                               base_dir::AbstractString = INT_DIR,
                               imp_dir::AbstractString  = INT_DIR,
                               overwrite::Bool          = false,
                               most_known_candidates::Vector{String} = String[])

    boots = load_all_bootstraps(; years, dir = base_dir, quiet = true)
    n     = length(boots)

    prog  = pm.Progress(n; desc = "Imputing bootstraps", barlen = 30)

    out   = OrderedCollections.OrderedDict{Int,String}()

    for (yr, bt) in boots
        @info "Imputing year $(yr)…"
        out[yr] = impute_and_save(bt;
                                  dir           = imp_dir,
                                  overwrite     = overwrite,
                                  most_known_candidates = most_known_candidates)
        pm.next!(prog)                  # advance progress bar
    end

    return out
end
 =#
#= 
function impute_all_bootstraps(; years    = nothing,
                               base_dir::AbstractString = INT_DIR,
                               imp_dir ::AbstractString = INT_DIR,
                               overwrite::Bool          = false,
                               most_known_candidates::Vector{String} = String[])

    # ── discover candidate files ------------------------------------------------
    rx = r"boot_(\d{4})\.jld2$"                # capture the 4-digit year
    files = readdir(base_dir; join = true)

    worklist = Vector{Tuple{Int,String}}()     # (year, path)

    for p in files
        m = match(rx, basename(p))
        m === nothing && continue
        yr = parse(Int, m.captures[1])

        if years !== nothing                    # filter if user requested
            yrs = isa(years, Integer) ? (years,) : years
            yr ∉ yrs && continue
        end

        push!(worklist, (yr, p))
    end

    isempty(worklist) && error("No bootstrap files found matching $rx in $base_dir")
    sort!(worklist; by = first)                 # chronological order

    # ── progress bar ------------------------------------------------------------
    prog = pm.Progress(length(worklist); desc = "Imputing bootstraps", barlen = 30)
    out  = OrderedDict{Int,String}()

    # ── streaming loop ----------------------------------------------------------
    for (yr, path) in worklist
        # load the two variables that save_bootstrap wrote
        reps = cfg = nothing
        @load path reps cfg
        bt = (data = reps, cfg = cfg, path = path)

        @info "Imputing year $yr…"

        out[yr] = impute_and_save(bt;
                                  dir               = imp_dir,
                                  overwrite         = overwrite,
                                  most_known_candidates = most_known_candidates)

        # cleanup
        bt = reps = cfg = nothing
        GC.gc()
        pm.next!(prog)
    end

    return out
end =#

struct ImputedYear
    year::Int
    # Dict(:zero => [path1, path2, …], :random => …, :mice => …)
    paths::Dict{Symbol,Vector{String}}
end

"""
    getrep(iy::ImputedYear, variant::Symbol, i::Int) -> DataFrame

Load the *i*-th replicate of `variant` for that year.
"""
function getrep(iy::ImputedYear, variant::Symbol, i::Int)
    p = iy.paths[variant][i]
    df = nothing; @load p df        # `df` is how we store it below
    return df
end

Base.getindex(iy::ImputedYear, variant::Symbol, i::Int) = getrep(iy, variant, i)

const IMP_DATA_DIR = joinpath(INT_DIR, "imputed_data"); mkpath(IMP_DATA_DIR)

function impute_bootstrap_to_files(path_boot::String;
                                   imp_dir::AbstractString = IMP_DATA_DIR,
                                   overwrite::Bool         = false,
                                   most_known_candidates   = String[])

    reps = cfg = nothing
    @load path_boot reps cfg                 # same vars saved by save_bootstrap
    year = cfg.year

    # ---------------- per-variant path collectors -----------------
    var_syms = (:zero, :random, :mice)
    paths_dict = Dict(var => Vector{String}(undef, length(reps)) for var in var_syms)

    for (i, df_raw) in enumerate(reps)
        imp = imputation_variants(df_raw, cfg.candidates, cfg.demographics;
                                  most_known_candidates)

        for var in var_syms
            file = joinpath(imp_dir,
                    "imp_$(year)_rep$(i)_$(String(var)).jld2")
            if !overwrite && isfile(file)
                @warn "reusing $(file)"
            else
                df = imp[var]                                # DataFrame
                @save file df
            end
            paths_dict[var][i] = file
        end

        # ---------- free memory for this replicate ----------
        imp = reps[i] = df_raw = nothing
        GC.gc()
    end

    # tiny index object
    index = ImputedYear(year, paths_dict)
    ind_file = joinpath(imp_dir, "index_$(year).jld2")
    @save ind_file index

    @info "Finished imputation for $(year) → $(ind_file)"
    return ind_file
end


function impute_all_bootstraps(; years = nothing,
                               base_dir = INT_DIR,
                               imp_dir  = IMP_DATA_DIR,
                               overwrite = false,
                               most_known_candidates = String[])

    rx = r"boot_(\d{4})\.jld2$"
    files = filter(p -> occursin(rx, basename(p)), readdir(base_dir; join=true))
    isempty(files) && error("No bootstrap files found in $(base_dir)")

    wanted = years === nothing ? nothing :
             isa(years,Integer) ? Set([years]) : Set(years)

    worklist = Tuple{Int,String}[]
    for p in files
        yr = parse(Int, match(rx, basename(p)).captures[1])
        (wanted !== nothing && yr ∉ wanted) && continue
        push!(worklist, (yr, p))
    end
    sort!(worklist; by = first)

    prog = pm.Progress(length(worklist); desc = "Imputing bootstraps", barlen = 30)
    out  = OrderedDict{Int,String}()

    for (yr, p) in worklist
        @info "Imputing year $yr …"
        out[yr] = impute_bootstrap_to_files(p;
                     imp_dir = imp_dir,
                     overwrite = overwrite,
                     most_known_candidates = most_known_candidates)
        GC.gc()
        pm.next!(prog)
    end
    return out
end


# ---------------------------------------------------------------------
# directory that will hold   imp_YYYY_repN_variant.jld2   and index_YYYY.jld2

# ---------------------------------------------------------------------
# 1 · impute ONE year's bootstrap and stream each replicate to its own file
# ---------------------------------------------------------------------
#= function _impute_year_to_files(reps::Vector{DataFrame},
                               cfg::ElectionConfig;
                               imp_dir::AbstractString = IMP_DATA_DIR,
                               overwrite::Bool = false,
                               most_known_candidates = String[])

    year  = cfg.year
    nboot = length(reps)
    variants = (:zero, :random, :mice)

    # Dict(:zero => Vector{String}(nboot), …)
    paths = Dict(var => Vector{String}(undef, nboot) for var in variants)

    for i in 1:nboot
        df_raw = reps[i]                       # original replicate
        imp    = imputation_variants(df_raw,
                                     cfg.candidates,
                                     cfg.demographics;
                                     most_known_candidates)

        for var in variants
            file = joinpath(imp_dir, "imp_$(year)_rep$(i)_$(String(var)).jld2")
            if overwrite || !isfile(file)
                df = imp[var]                  # DataFrame for that variant
                @save file df
            end
            paths[var][i] = file
        end

        # ---- drop local references; leave `reps` untouched -------------
        imp = df_raw = nothing
        GC.gc()
    end

    # tiny index object
    index   = ImputedYear(year, paths)
    idxfile = joinpath(imp_dir, "index_$(year).jld2")
    @save idxfile index
    return idxfile
end =#

function _impute_year_to_files(reps::Vector{DataFrame},
                               cfg::ElectionConfig;
                               imp_dir::AbstractString = IMP_DATA_DIR,
                               overwrite::Bool = false,
                               most_known_candidates = String[])

    year      = cfg.year
    idxfile   = joinpath(imp_dir, "index_$(year).jld2")
    variants  = (:zero, :random, :mice)
    nboot     = length(reps)

    # ────────────────────────────────────────────────────────────────────
    # 1. Fast path: cached index exists
    # ────────────────────────────────────────────────────────────────────
    if !overwrite && isfile(idxfile)
        @info "Reusing cached imputed index for year $year → $(idxfile)"
        return idxfile
    end

    # helper that lists all expected replicate files
    expected_path(var, i) = joinpath(
        imp_dir, "imp_$(year)_rep$(i)_$(String(var)).jld2")

    # ────────────────────────────────────────────────────────────────────
    # 2.  Could we build the index without recomputation?
    # ────────────────────────────────────────────────────────────────────
    if !overwrite
        all_exist = true
        for i in 1:nboot, var in variants
            isfile(expected_path(var, i)) || (all_exist = false; break)
        end
        if all_exist
            paths = Dict(var => [expected_path(var, i) for i in 1:nboot]
                         for var in variants)
            index = ImputedYear(year, paths)
            @save idxfile index
            @info "Rebuilt index for year $year without re-imputation."
            return idxfile
        end
    end

    # ────────────────────────────────────────────────────────────────────
    # 3.  Run full imputation (some files missing or overwrite=true)
    # ────────────────────────────────────────────────────────────────────
    @info "Running imputation for year $year …"
    paths = Dict(var => Vector{String}(undef, nboot) for var in variants)

    for i in 1:nboot
        df_raw = reps[i]
        imp    = imputation_variants(df_raw,
                                     cfg.candidates,
                                     cfg.demographics;
                                     most_known_candidates)

        for var in variants
            file = expected_path(var, i)
            if overwrite || !isfile(file)
                df = imp[var]
                @save file df
            end
            paths[var][i] = file
        end

        imp    = df_raw = nothing        # local cleanup
        GC.gc()
    end

    index = ImputedYear(year, paths)
    @save idxfile index
    @info "Saved imputed index for year $year → $(idxfile)"
    return idxfile
end


# ---------------------------------------------------------------------
# 2 · top-level driver starting from your in-memory `f3`
# ---------------------------------------------------------------------
function impute_from_f3(f3::OrderedDict;
                        years = nothing,
                        imp_dir::AbstractString = IMP_DATA_DIR,
                        overwrite::Bool = false,
                        most_known_candidates = String[])

    wanted = years === nothing        ? sort(collect(keys(f3))) :
             isa(years,Integer)       ? [years]                 :
             sort(collect(years))

    prog = pm.Progress(length(wanted); desc = "Imputing bootstraps", barlen = 30)
    out  = OrderedDict{Int,String}()

    for yr in wanted
        entry = f3[yr]                       # (data = reps, cfg = cfg, path = …)
        reps  = entry.data
        cfg   = entry.cfg

        @info "Imputing year $yr …"
        out[yr] = _impute_year_to_files(reps, cfg;
                                        imp_dir    = imp_dir,
                                        overwrite  = overwrite,
                                        most_known_candidates = most_known_candidates)

        GC.gc()                             # reclaim before next year
        pm.next!(prog)
    end
    return out
end

const IMP_PREFIX  = "boot_imp_"   # change here if you rename files
const IMP_DIR     = INT_DIR       # default directory to look in

# ————————————————————————————————————————————————————————————————
# 1.  Single-year loader
# ————————————————————————————————————————————————————————————————
"""
    load_imputed_bootstrap(year;
                           dir   = IMP_DIR,
                           quiet = false)  -> NamedTuple

Load `dir/boot_imp_YEAR.jld2` and return the stored NamedTuple
`(data = Dict, cfg = ElectionConfig, path = String)`.
"""
function load_imputed_bootstrap(year::Integer;
                                dir::AbstractString = IMP_DIR,
                                quiet::Bool = false)

    path = joinpath(dir, "$(IMP_PREFIX)$(year).jld2")
    isfile(path) || error("File not found: $(path)")

    imp = nothing
    @load path imp
    !quiet && @info "Loaded imputed bootstrap for year $(year) ← $(path)"
    return imp
end

function load_all_imputed_bootstraps(; years  = nothing,
                                     dir::AbstractString = IMP_DIR,
                                     quiet::Bool = false)

    # 1 — discover candidate files ------------------------------------------------
    allfiles = readdir(dir; join = true)
    paths = filter(p -> startswith(basename(p), IMP_PREFIX) &&
                    endswith(p, ".jld2"), allfiles)

    isempty(paths) && error("No imputed bootstrap files found in $(dir)")

    # 2 — decide which years we actually want ------------------------------------
    wanted = years === nothing       ? nothing :
             isa(years,Integer)      ? Set([years]) :
             Set(years)

    # 3 — load, build OrderedDict -------------------------------------------------
    out = OrderedCollections.OrderedDict{Int,NamedTuple}()

    for p in sort(paths)                           # alphabetical == chronological
        yr = parse(Int, splitext(basename(p)[length(IMP_PREFIX)+1:end])[1])
        (wanted !== nothing && !(yr in wanted)) && continue

        imp = nothing
        @load p imp
        !quiet && @info "Loaded imputed bootstrap $(yr) ← $(p)"
        out[yr] = imp
    end

    return out
end


function load_imputed_year(year::Int;
                           dir::AbstractString = IMP_DATA_DIR)::ImputedYear
    idxfile = joinpath(dir, "index_$(year).jld2")
    isfile(idxfile) || error("index file not found: $(idxfile)")
    return JLD2.load(idxfile, "index")    # returns ImputedYear struct
end

# TODO: later, write a variant that takes just imp and config
# and loads f3 from disk, cakculate the sets, cleans it from disk, and proceeds


const CANDLOG = joinpath(INT_DIR, "candidate_set_warnings.log")
#= 
function generate_profiles_for_year(year::Int,
                                    f3_entry::NamedTuple,
                                    imps_entry::NamedTuple)

    cfg            = f3_entry.cfg
    reps_raw       = f3_entry.data             # Vector{DataFrame}
    variants_dict  = imps_entry.data           # Dict{Symbol,Vector{DataFrame}}
    m_values       = cfg.m_values_range

    result = OrderedDict{String, OrderedDict{Int,OrderedDict{Symbol,Vector{DataFrame}}}}()

    for scen in cfg.scenarios
        # —————————————————————————————————————————————————————
        # 1) compute the “full” candidate set once (m = max_candidates)
        sets = unique(map(df ->
            compute_candidate_set(df;
                candidate_cols = cfg.candidates,
                m              = cfg.max_candidates,
                force_include  = scen.candidates),
          reps_raw))

        if length(sets) != 1
            msg = "Year $year, scenario $(scen.name): " *
                  "found $(length(sets)) distinct candidate sets; " *
                  "using the first. Sets: $sets"
            @warn msg

            # append to disk log
            open(CANDLOG, "a") do io
                ts = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
                println(io, "[$ts] $msg")
            end
        end

        full_list = sets[1]   # pick the first, regardless

        # —————————————————————————————————————————————————————
        # 2) for each m, trim & build profile DFs
        m_map = OrderedDict{Int, OrderedDict{Symbol,Vector{DataFrame}}}()
        for m in m_values
            trimmed = first(full_list, m)
            var_map = OrderedDict{Symbol,Vector{DataFrame}}()
            for (variant, reps_imp) in variants_dict
                profiles = Vector{DataFrame}(undef, length(reps_imp))
                for (i, df_imp) in enumerate(reps_imp)
                    profiles[i] = profile_dataframe(
                                      df_imp;
                                      score_cols = trimmed,
                                      demo_cols  = cfg.demographics)
                end
                var_map[variant] = profiles
            end
            m_map[m] = var_map
        end

        result[scen.name] = m_map
    end

    return result
end
 =#



function generate_profiles_for_year(year::Int,
                                    f3_entry::NamedTuple,
                                    imps_entry::NamedTuple)

    cfg            = f3_entry.cfg
    reps_raw       = f3_entry.data
    variants_dict  = imps_entry.data
    m_values       = cfg.m_values_range

    result = OrderedDict{String,OrderedDict{Int,OrderedDict{Symbol,Vector{DataFrame}}}}()

    for scen in cfg.scenarios
        sets = unique(map(df ->
            compute_candidate_set(df;
                candidate_cols = cfg.candidates,
                m              = cfg.max_candidates,
                force_include  = scen.candidates),
            reps_raw))

        length(sets) != 1 && @warn "Year $year, scenario $(scen.name): $(length(sets)) candidate sets; using first."
        full_list = sets[1]

        m_map = OrderedDict{Int,OrderedDict{Symbol,Vector{DataFrame}}}()

        for m in m_values
            trimmed  = Symbol.(first(full_list, m))          # ordered Vector{String}
            var_map  = OrderedDict{Symbol,Vector{DataFrame}}()

            for (variant, reps_imp) in variants_dict
                profiles = Vector{DataFrame}(undef, length(reps_imp))

                for (i, df_imp) in enumerate(reps_imp)
                    df = profile_dataframe(
                             df_imp;
                             score_cols = trimmed,
                             demo_cols  = cfg.demographics)
                    compress_rank_column!(df, trimmed; col = :profile)
                    # ------------- NEW LINE ----------------------------------
                    metadata!(df, "candidates", Symbol.(trimmed))  
                    # ---------------------------------------------------------

                    profiles[i] = df
                end
                var_map[variant] = profiles
            end
            m_map[m] = var_map
        end
        result[scen.name] = m_map
    end
    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# directories / tiny types (unchanged)
# ─────────────────────────────────────────────────────────────────────────────
const PROFILES_DATA_DIR = joinpath(INT_DIR, "profiles_data")
mkpath(PROFILES_DATA_DIR)

struct ProfilesSlice
    year::Int
    scenario::String
    m::Int
    cand_list::Vector{Symbol}              # ordered for (en/de)coding
    paths::Dict{Symbol,Vector{String}}     # variant ⇒ file paths
end

Base.getindex(ps::ProfilesSlice, var::Symbol, i::Int) = begin
    p = ps.paths[var][i]; df = nothing; JLD2.@load p df; df
end


function generate_profiles_for_year_streamed_from_index(
            year::Int,
            f3_entry::NamedTuple,
            iy::ImputedYear;
            out_dir::AbstractString = PROFILES_DATA_DIR,
            overwrite::Bool         = false)

    cfg            = f3_entry.cfg
    reps_raw       = f3_entry.data
    m_values       = cfg.m_values_range
    variants       = collect(keys(iy.paths))           # e.g. (:zero,:random,:mice)
    n_by_var       = Dict(v => length(iy.paths[v]) for v in variants)

    result = OrderedDict{String,OrderedDict{Int,ProfilesSlice}}()

    for scen in cfg.scenarios
        sets = unique(map(df ->
            compute_candidate_set(df;
                candidate_cols = cfg.candidates,
                m              = cfg.max_candidates,
                force_include  = scen.candidates),
            reps_raw))
        length(sets) != 1 && @warn "Year $year, scenario $(scen.name): " *
                                   "$(length(sets)) candidate sets; using the first."
        full_cset = sets[1]

        scen_map = OrderedDict{Int,ProfilesSlice}()

        for m in m_values
            cand_syms  = Symbol.(first(full_cset, m))
            paths_prof = Dict(v => Vector{String}(undef, n_by_var[v]) for v in variants)

            rep_counter = 0      # throttle GC

            for var in variants
                n_rep = n_by_var[var]

                for i in 1:n_rep
                    fprof = joinpath(out_dir,
                             "prof_$(year)_$(scen.name)_m$(m)_rep$(i)_" *
                             "$(String(var)).jld2")

                    # -------- fast‑skip if file already present --------------
                    if !overwrite && isfile(fprof)
                        paths_prof[var][i] = fprof
                        @debug "exists, skipping $(basename(fprof))"
                        continue
                    end

                    # -------- otherwise build & save -------------------------
                    df_imp = iy[var, i]

                    df = profile_dataframe(df_imp;
                            score_cols = cand_syms,
                            demo_cols  = cfg.demographics)
                    compress_rank_column!(df, cand_syms; col = :profile)
                    metadata!(df, "candidates", cand_syms)

                    JLD2.@save fprof df
                    @info "writing $(basename(fprof))"
                    paths_prof[var][i] = fprof

                    df = df_imp = nothing
                    rep_counter += 1
                    rep_counter % 10 == 0 && GC.gc()
                end
            end

            slice = ProfilesSlice(year, scen.name, m, cand_syms, paths_prof)
            scen_map[m] = slice
        end
        result[scen.name] = scen_map
    end

    idxfile = joinpath(out_dir, "profiles_index_$(year).jld2")
    JLD2.@save idxfile result
    @info "Encoded profiles for $year written; index at $(idxfile)"
    return result
end




#= 
"""Lazy loader that accepts a DataFrame *or* a reference."""
load_df(ref) = ref isa DataFrame ? ref : load_dataframe(ref)  # write this once

function generate_profiles_for_year(year,
                                    f3_entry,
                                    imps_entry)

    cfg            = f3_entry.cfg
    reps_raw_refs  = f3_entry.data             # Vector{DataFrame|Handle}
    variants_refs  = imps_entry.data           # Dict{Symbol,Vector{DataFrame|Handle}}
    m_values       = cfg.m_values_range

    result = OrderedDict{String, OrderedDict{Int,OrderedDict{Symbol,Vector{DataFrame}}}}()

    for scen in cfg.scenarios
        # ────────────────────────────────────────────────────────────────
        # 1) Determine (once) the reference candidate set, streaming rep-by-rep
        # ────────────────────────────────────────────────────────────────
        ref_set   = nothing
        divergent = false

        for ref in reps_raw_refs
            df   = load_df(ref)                                   # ← lazy fetch
            cset = compute_candidate_set(df;
                       candidate_cols = cfg.candidates,
                       m              = cfg.max_candidates,
                       force_include  = scen.candidates)

            if ref_set === nothing
                ref_set = cset
            elseif ref_set != cset
                divergent = true
                break
            end
            GC.gc()                                               # drop df ASAP
        end

        if divergent
            sets = map(refs -> compute_candidate_set(load_df(refs);
                                  candidate_cols = cfg.candidates,
                                  m              = cfg.max_candidates,
                                  force_include  = scen.candidates),
                       first(reps_raw_refs, 7))                    # cheap summary
            msg = "Year $year, scenario $(scen.name): found ≥2 distinct candidate sets; using the first. Sets: $(sets)"
            @warn msg
            open(CANDLOG, "a") do io
                println(io, "[$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))] $msg")
            end
        end

        full_list = ref_set
        # ────────────────────────────────────────────────────────────────
        # 2) Build the profiles one replicate at a time
        # ────────────────────────────────────────────────────────────────
        m_map = OrderedDict{Int, OrderedDict{Symbol,Vector{DataFrame}}}()

        for m in m_values
            trimmed   = first(full_list, m)
            var_map   = OrderedDict{Symbol,Vector{DataFrame}}()

            for (variant, reps_imp_refs) in variants_refs
                profiles = Vector{DataFrame}(undef, length(reps_imp_refs))

                for (i, ref_imp) in enumerate(reps_imp_refs)
                    df_imp        = load_df(ref_imp)              # ← lazy fetch
                    profiles[i]   = profile_dataframe(
                                         df_imp;
                                         score_cols = trimmed,
                                         demo_cols  = cfg.demographics)
                    GC.gc()                                       # drop df_imp
                end
                var_map[variant] = profiles
            end

            m_map[m] = var_map
            GC.gc()                                               # reclaim temps
        end

        result[scen.name] = m_map
        GC.gc()
    end

    return result
end

 =#

function generate_all_profiles(f3::OrderedDict{Int,NamedTuple},
                               imps::OrderedDict{Int,NamedTuple};
                               years = nothing)

    wanted = years === nothing       ? keys(f3) :
             isa(years, Integer)      ? [years]      :
             years

    all_profiles = OrderedDict{Int,Any}()
    for yr in sort(collect(wanted))
        f3_entry  = f3[yr]
        imps_entry = imps[yr]
        @info "Building profiles for year $yr"
        all_profiles[yr] = generate_profiles_for_year(yr, f3_entry, imps_entry)
    end
    return all_profiles
end



const PROFILE_FILE = joinpath(INT_DIR, "all_profiles.jld2")

"""
    save_or_load_all_profiles(f3, imps;
                              path         = PROFILE_FILE,
                              overwrite    = false,
                              verbose      = true,
                              years        = nothing)

If `path` already exists *and* `overwrite=false`, emits a warning,
loads the saved `profiles` object, and returns it.

Otherwise, calls `generate_all_profiles(f3, imps; years=years)`,
saves the result to `path`, and returns it.
"""
function save_or_load_all_profiles(f3, imps;
                                   path      ::AbstractString = PROFILE_FILE,
                                   overwrite ::Bool           = false,
                                   verbose   ::Bool           = true,
                                   years               = nothing)

    # ensure directory
    mkpath(dirname(path))

    if isfile(path) && !overwrite
        verbose && @warn "Profiles file already exists at $path; loading it instead of regenerating."
        profiles = nothing
        @load path profiles
        return profiles
    end

    # generate + save
    verbose && @info "Generating all profiles (this may take a while)…"
    profiles = generate_all_profiles(f3, imps; years = years)

    @save path profiles
    verbose && @info "Saved all profiles to $path"
    return profiles
end


"""
    load_all_profiles(; path = PROFILE_FILE, verbose = true) -> OrderedDict

Load the `profiles` object from disk at `path`. Throws an error if the file
does not exist.
"""
function load_all_profiles(; path    ::AbstractString = PROFILE_FILE,
                            verbose ::Bool           = true)

    isfile(path) || error("No profiles file found at $path -- run `save_or_load_all_profiles` first.")
    profiles = nothing
    @load path profiles
    verbose && @info "Loaded profiles from $path"
    return profiles
end




const PROFILE_DIR = joinpath(INT_DIR, "profiles")
mkpath(PROFILE_DIR)   # ensure it exists

# ────────────────────────────────────────────────────────────────────────────────
"""
    save_or_load_profiles_for_year(year, f3, imps;
                                  dir       = PROFILE_DIR,
                                  overwrite = false,
                                  verbose   = true)

For the given `year`, if `dir/profiles_YEAR.jld2` exists and `overwrite=false`,
issues a warning and loads it.  Otherwise:

  • Calls `generate_profiles_for_year(year, f3[year], imps[year])`  
  • Saves the result as `profiles` in `dir/profiles_YEAR.jld2`  
  • Returns the `profiles` object.
"""
function save_or_load_profiles_for_year(year::Int,
                                        f3,
                                        imps;
                                        dir::AbstractString = PROFILE_DIR,
                                        overwrite::Bool     = false,
                                        verbose::Bool       = true)

    path = joinpath(dir, "profiles_$(year).jld2")

    if isfile(path) && !overwrite
        verbose && @warn "Profiles for $year already exist at $path; loading cache."
        profiles = nothing
        @load path profiles
        return profiles
    end

    verbose && @info "Generating profiles for year $year…"
    profiles = generate_profiles_for_year(year, f3[year], imps[year])
    @save path profiles
    verbose && @info "Saved profiles for year $year → $path"
    return profiles
end

"""
    save_or_load_all_profiles_per_year(f3, imps;
                                       years     = nothing,
                                       dir::AbstractString = PROFILE_DIR,
                                       overwrite ::Bool     = false,
                                       verbose   ::Bool     = true)
                                       
Iterate over each `year` in `f3` (or the subset `years`) and call
`save_or_load_profiles_for_year`.  Returns an `OrderedDict{Int,profiles}`.
"""
function save_or_load_all_profiles_per_year(f3,
                                            imps;
                                            years     = nothing,
                                            dir::AbstractString = PROFILE_DIR,
                                            overwrite ::Bool     = false,
                                            verbose   ::Bool     = true)

    wanted = years === nothing       ? sort(collect(keys(f3))) :
             isa(years, Integer)      ? [years]              :
             sort(collect(years))

    all_profiles = OrderedDict{Int,Any}()
    for yr in wanted
        if !haskey(f3, yr)
            @warn "No bootstrap for year $yr; skipping"
            continue
        end
        all_profiles[yr] = save_or_load_profiles_for_year(
                               yr, f3, imps;
                               dir       = dir,
                               overwrite = overwrite,
                               verbose   = verbose)
    end

    return all_profiles
end

# ────────────────────────────────────────────────────────────────────────────────
"""
    load_profiles_for_year(year;
                            dir     = PROFILE_DIR,
                            verbose = true) -> profiles

Loads `dir/profiles_YEAR.jld2` and returns the stored `profiles` object.
Errors if missing.
"""
function load_profiles_for_year(year::Int;
                                dir::AbstractString = PROFILE_DIR,
                                verbose::Bool       = true)

    path = joinpath(dir, "profiles_$(year).jld2")
    isfile(path) || error("No profiles JLD2 found for $year at $path")
    profiles = nothing
    @load path profiles
    verbose && @info "Loaded profiles for year $year ← $path"
    return profiles
end

"Make a Dict<measure,<variant,Vector>> skeleton with empty vectors."
function init_accumulator(var_syms, measure_syms)
    accum = Dict{Symbol,Dict{Symbol,Vector{Float64}}}()
    for meas in measure_syms
        inner = Dict{Symbol,Vector{Float64}}()
        for var in var_syms
            inner[var] = Float64[]          # will push! into it
        end
        accum[meas] = inner
    end
    return accum
end

"Append values from `meas_one_rep` (1‑replicate output) into `accum`."
@inline function update_accumulator!(accum, meas_one_rep)
    for (meas, vdict) in meas_one_rep           # meas ⇒ variant ⇒ Vector(1)
        inner = accum[meas]
        for (var, vec1) in vdict
            push!(inner[var], vec1[1])          # vec1 has length 1
        end
    end
    return
end

#= function apply_measures_for_year(
    profiles_year::OrderedDict{String,<:Any}
)::OrderedDict{String,OrderedDict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}}

    out = OrderedDict{String,OrderedDict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}}()

    for (scen, m_map) in profiles_year
        scen_out = OrderedDict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}()

        for (m, slice) in m_map               # slice :: ProfilesSlice
            @assert slice isa ProfilesSlice   # streaming logic expects this

            variants  = collect(keys(slice.paths))
            n_rep_max = maximum(length(slice.paths[v]) for v in variants)

            # ── first replicate: discover measure names ─────────────
            var_map1 = Dict(var => [slice[var,1]] for var in variants)
            decode_each!(var_map1)
            meas1    = apply_all_measures_to_bts(var_map1)
            meas_syms = collect(keys(meas1))

            accum = init_accumulator(variants, meas_syms)
            update_accumulator!(accum, meas1)

            # ── remaining replicates ────────────────────────────────
            rep_counter = 1
            for i in 2:n_rep_max
                var_map = Dict{Symbol,Vector{DataFrame}}()

                for var in variants
                    length(slice.paths[var]) < i && continue
                    df = slice[var, i]
                    decode_profile_column!(df)
                    var_map[var] = [df]
                end
                isempty(var_map) && continue

                meas_i = apply_all_measures_to_bts(var_map)
                update_accumulator!(accum, meas_i)

                rep_counter += 1
                rep_counter % 10 == 0 && GC.gc()
            end

            scen_out[m] = accum
            GC.gc()
        end
        out[scen] = scen_out
    end
    return out
end
 =#


function apply_measures_for_year(
    profiles_year::OrderedDict{String,<:Any}
)::OrderedDict{String,OrderedDict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}}

    out = OrderedDict{String,OrderedDict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}}()

    for (scen, m_map) in profiles_year
        scen_out = OrderedDict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}()

        for (m, slice) in m_map            # slice :: ProfilesSlice
            @assert slice isa ProfilesSlice

            variants   = collect(keys(slice.paths))
            n_rep_max  = maximum(length(slice.paths[v]) for v in variants)

            #####  first replicate: discover measure names  #####
            var_map1 = Dict(var => [slice[var, 1]] for var in variants)
            decode_each!(var_map1)
            meas1     = apply_all_measures_to_bts(var_map1)
            meas_syms = collect(keys(meas1))

            accum = init_accumulator(variants, meas_syms)
            update_accumulator!(accum, meas1)

            #####  progress bar  #####
            prog = pm.Progress(n_rep_max - 1;  desc = "[$scen|m=$m]", barlen = 30)

            #####  remaining replicates  #####
            rep_counter = 1
            for i in 2:n_rep_max
                var_map = Dict{Symbol,Vector{DataFrame}}()

                for var in variants
                    length(slice.paths[var]) < i && continue
                    df = slice[var, i]
                    decode_profile_column!(df)
                    var_map[var] = [df]
                end
                isempty(var_map) && continue

                meas_i = apply_all_measures_to_bts(var_map)
                update_accumulator!(accum, meas_i)

                pm.next!(prog)                           # advance bar
                rep_counter += 1
                rep_counter % 10 == 0 && GC.gc()
            end
            pm.finish!(prog)

            scen_out[m] = accum
            GC.gc()
        end
        out[scen] = scen_out
    end
    return out
end

#= 
"""
    apply_measures_for_year(profiles_year::OrderedDict)

Given `profiles_year` produced by `save_or_load_profiles_for_year`, return

    OrderedDict{String,OrderedDict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}}

so that

    result[scenario][m][measure][variant]  # → Vector{Float64}

Internally it just does:

    measure_map = apply_all_measures_to_bts(var_map)

for each (scenario → m → var_map).
"""
#= function apply_measures_for_year(
    profiles_year::OrderedDict{String,<:Any}
)::OrderedDict{String,OrderedDict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}}

    out = OrderedDict{String,OrderedDict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}}()

    for (scen, m_map) in profiles_year
        scen_out = OrderedDict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}()
        for (m, var_map) in m_map
            # var_map :: Dict{Symbol, Vector{DataFrame}}
            measures = apply_all_measures_to_bts(var_map)
            scen_out[m] = measures
        end
        out[scen] = scen_out
    end

    return out
end
 =# =#


function apply_measures_all_years(
    profiles_all::Dict{Int,Any};
    years = nothing
)::OrderedDict{Int,OrderedDict{String,OrderedDict{Int,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}}}

    wanted = years === nothing       ? sort(collect(keys(profiles_all))) :
             isa(years, Integer)      ? [years]                      :
             sort(collect(years))

    all_out = OrderedDict{Int,Any}()
    for yr in wanted
        haskey(profiles_all, yr) || continue
        @info "Applying measures for year $yr"
        all_out[yr] = apply_measures_for_year(profiles_all[yr])
    end
    return all_out
end


const GLOBAL_MEASURE_DIR = joinpath(INT_DIR, "global_measures")
mkpath(GLOBAL_MEASURE_DIR)   # ensure the directory exists

"""
    save_or_load_measures_for_year(year, profiles_year;
                                   dir       = GLOBAL_MEASURE_DIR,
                                   overwrite = false,
                                   verbose   = true)

For a single `year`:

- If `dir/measures_YEAR.jld2` exists and `overwrite == false`, emits a warning and loads `measures` from disk.
- Otherwise, runs `apply_measures_for_year(profiles_year)`, saves the result under the name `measures`, and returns it.
"""
function save_or_load_measures_for_year(year,
                                        profiles_year;
                                        dir::AbstractString = GLOBAL_MEASURE_DIR,
                                        overwrite::Bool     = false,
                                        verbose::Bool       = true)

    path = joinpath(dir, "measures_$(year).jld2")

    if isfile(path) && !overwrite
        verbose && @warn "Global measures for $year already cached at $path; loading."
        measures = nothing
        @load path measures
        return measures
    end

    verbose && @info "Computing global measures for year $year…"
    measures = apply_measures_for_year(profiles_year)
    @save path measures
    verbose && @info "Saved global measures for year $year → $path"
    return measures
end

"""
    save_or_load_all_measures_per_year(profiles_all;
                                      years     = nothing,
                                      dir       = GLOBAL_MEASURE_DIR,
                                      overwrite = false,
                                      verbose   = true)

Loop over each `year` in `profiles_all` (or the subset `years`) and call
`save_or_load_measures_for_year`. Returns an `OrderedDict{Int,measures}`.
"""
function save_or_load_all_measures_per_year(profiles_all;
                                            years     = nothing,
                                            dir::AbstractString = GLOBAL_MEASURE_DIR,
                                            overwrite ::Bool     = false,
                                            verbose   ::Bool     = true)

    wanted = years === nothing       ? sort(collect(keys(profiles_all))) :
             isa(years, Integer)      ? [years]              :
             sort(collect(years))

    all_measures = OrderedDict{Int,Any}()
    for yr in wanted
        if !haskey(profiles_all, yr)
            @warn "No profiles for year $yr; skipping measures."
            continue
        end
        all_measures[yr] = save_or_load_measures_for_year(
                               yr,
                               profiles_all[yr];
                               dir       = dir,
                               overwrite = overwrite,
                               verbose   = verbose)
    end

    return all_measures
end

"""
    load_measures_for_year(year;
                          dir     = GLOBAL_MEASURE_DIR,
                          verbose = true) -> measures

Load `dir/measures_YEAR.jld2` and return the `measures` object.
Errors if missing.
"""
function load_measures_for_year(year::Int;
                                dir::AbstractString = GLOBAL_MEASURE_DIR,
                                verbose::Bool       = true)

    path = joinpath(dir, "measures_$(year).jld2")
    isfile(path) || error("No global measures file for $year at $path")
    measures = nothing
    @load path measures
    verbose && @info "Loaded global measures for year $year ← $path"
    return measures
end 


#= 


function apply_group_metrics_for_year(profiles_year, cfg)
    # profiles_year: any dict-like scen ⇒ m ⇒ variant ⇒ Vector{DF}
    # cfg.demographics: Vector of String or Symbol
    out = OrderedDict()

    for (scen, m_map) in profiles_year
        scen_out = OrderedDict()
        for (m, var_map) in m_map
            dem_out = OrderedDict()
            for dem in cfg.demographics
                # ensure we pass a Symbol to your metrics function
                dem_sym = dem isa Symbol ? dem : Symbol(dem)
                @info "  → scenario=$scen, m=$m, demographic=$dem_sym"
                dem_out[dem_sym] = bootstrap_group_metrics(var_map, dem_sym)
            end
            scen_out[m] = dem_out
        end
        out[scen] = scen_out
    end

    return out
end


function apply_group_metrics_all_years(profiles_all, f3; years=nothing)
    # profiles_all: year ⇒ per‐year profiles (scenario ⇒ m ⇒ variant ⇒ DF)
    # f3:          year ⇒ (data,cfg,path)
    wanted = years === nothing ? sort(collect(keys(profiles_all))) :
             isa(years, Integer) ? [years] :
             sort(collect(years))

    all_out = OrderedDict()
    for yr in wanted
        haskey(profiles_all, yr) || continue
        haskey(f3, yr)            || continue

        @info "Applying group metrics for year $yr"
        profiles_year = profiles_all[yr]
        cfg = f3[yr].cfg

        all_out[yr] = apply_group_metrics_for_year(profiles_year, cfg)
    end

    return all_out
end


const GROUP_DIR = joinpath(INT_DIR, "group_metrics")
mkpath(GROUP_DIR)   # ensure it exists

"""
    save_or_load_group_metrics_for_year(year,
                                        profiles_year,
                                        f3_entry;
                                        dir       = GROUP_DIR,
                                        overwrite = false,
                                        verbose   = true)

For a single `year`:

- If `dir/group_metrics_YEAR.jld2` exists and `overwrite==false`, emits a warning and loads `metrics` from disk.
- Otherwise, runs `apply_group_metrics_for_year(profiles_year, f3_entry.cfg)`, saves the result under the name `metrics`, and returns it.
"""
function save_or_load_group_metrics_for_year(year::Int,
                                             profiles_year,
                                             f3_entry;
                                             dir::AbstractString = GROUP_DIR,
                                             overwrite::Bool     = false,
                                             verbose::Bool       = true)

    path = joinpath(dir, "group_metrics_$(year).jld2")

    if isfile(path) && !overwrite
        verbose && @warn "Group metrics for $year already exist at $path; loading cache."
        metrics = nothing
        @load path metrics
        return metrics
    end

    verbose && @info "Computing group metrics for year $year…"
    metrics = apply_group_metrics_for_year(profiles_year, f3_entry.cfg)
    @save path metrics
    verbose && @info "Saved group metrics for year $year → $path"
    return metrics
end


"""
    save_or_load_all_group_metrics_per_year(profiles_all, f3;
                                            years     = nothing,
                                            dir       = GROUP_DIR,
                                            overwrite = false,
                                            verbose   = true)
                                            
Iterate over each `year` in `profiles_all` (or the subset `years`) and call
`save_or_load_group_metrics_for_year`. Returns an `OrderedDict{Int,metrics}`.
"""
function save_or_load_all_group_metrics_per_year(profiles_all,
                                                 f3;
                                                 years     = nothing,
                                                 dir::AbstractString = GROUP_DIR,
                                                 overwrite ::Bool     = false,
                                                 verbose   ::Bool     = true)

    wanted = years === nothing       ? sort(collect(keys(profiles_all))) :
             isa(years, Integer)      ? [years]              :
             sort(collect(years))

    all_metrics = OrderedDict{Int,Any}()
    for yr in wanted
        if !haskey(profiles_all, yr)
            @warn "No profiles for year $yr; skipping"
            continue
        end
        if !haskey(f3, yr)
            @warn "No bootstrap/config for year $yr; skipping"
            continue
        end

        all_metrics[yr] = save_or_load_group_metrics_for_year(
                              yr,
                              profiles_all[yr],
                              f3[yr];
                              dir       = dir,
                              overwrite = overwrite,
                              verbose   = verbose)
    end

    return all_metrics
end


"""
    load_group_metrics_for_year(year;
                                dir     = GROUP_DIR,
                                verbose = true) -> metrics

Load `dir/group_metrics_YEAR.jld2` and return the `metrics` object.
Errors if missing.
"""
function load_group_metrics_for_year(year::Int;
                                     dir::AbstractString = GROUP_DIR,
                                     verbose::Bool       = true)

    path = joinpath(dir, "group_metrics_$(year).jld2")
    isfile(path) || error("No group metrics file for $year at $path")
    metrics = nothing
    @load path metrics
    verbose && @info "Loaded group metrics for year $year ← $path"
    return metrics
end
 =#

const GROUP_DIR = joinpath(INT_DIR, "group_metrics"); mkpath(GROUP_DIR)

init_accum(vars, met_syms) =
    Dict(met => Dict(var => Float64[] for var in vars) for met in met_syms)



function update_accum!(accum::Dict, res::Dict, variants)
    for (met, vdict) in res
        inner = get!(accum, met) do
            # first time we see this metric → create inner dict with empty vectors
            Dict(var => Float64[] for var in variants)
        end
        for (var, vec1) in vdict          # vec1 length == 1
            push!(get!(inner, var, Float64[]), vec1[1])   # create variant slot if absent
        end
    end
end

# ─────────────────── streaming apply_group_metrics_for_year ───────────────────
function apply_group_metrics_for_year_streaming(
        profiles_year::OrderedDict{String,<:Any},
        cfg)

    out = OrderedDict()

    for (scen, m_map) in profiles_year
        scen_out = OrderedDict()

        for (m, slice) in m_map             # slice :: ProfilesSlice
            variants   = collect(keys(slice.paths))
            n_rep_max  = maximum(length(slice.paths[v]) for v in variants)

            dem_out = OrderedDict()

            for dem in cfg.demographics
                dem_sym = dem isa Symbol ? dem : Symbol(dem)
                @info "  → scenario=$scen, m=$m, dem=$dem_sym"

                accum = Dict{Symbol,Dict{Symbol,Vector{Float64}}}()
                prog  = pm.Progress(n_rep_max; desc="[$scen|m=$m|$dem_sym]", barlen=28)

                for i in 1:n_rep_max
                    var_map = Dict{Symbol,Vector{DataFrame}}()
                    for var in variants
                        length(slice.paths[var]) < i && continue
                        df = slice[var, i]
                        decode_profile_column!(df)
                        var_map[var] = [df]
                    end
                    isempty(var_map) && (next!(prog); continue)

                    res = bootstrap_group_metrics(var_map, dem_sym)
                    update_accum!(accum, res, variants)
                    pm.next!(prog)
                end
                pm.finish!(prog)
                dem_out[dem_sym] = accum
                GC.gc()
            end
            scen_out[m] = dem_out
        end
        out[scen] = scen_out
    end
    return out          # scenario ⇒ m ⇒ dem ⇒ metric ⇒ variant ⇒ Vector
end





#= 
function save_or_load_group_metrics_for_year(year::Int,
                                             profiles_year,
                                             f3_entry;
                                             dir::AbstractString = GROUP_DIR,
                                             overwrite::Bool     = false,
                                             verbose::Bool       = true)

    path = joinpath(dir, "group_metrics_$(year).jld2")

    if isfile(path) && !overwrite
        verbose && @warn "Group metrics for $year already cached at $path; loading."
        metrics = nothing; @load path metrics; return metrics
    end

    verbose && @info "Computing group metrics for year $year…"
    metrics = apply_group_metrics_for_year_streaming(profiles_year, f3_entry.cfg)
    @save path metrics
    verbose && @info "Saved group metrics for year $year → $path"
    return metrics
end


=#

 


# directory layout for per‑DataFrame caches
_perdf(dir, year, scen, m, dem, rep) =
    joinpath(dir, "per_df", string(year), string(scen), "m$m", string(dem),
             "rep$rep.jld2")

# ensure that a path’s parent directories exist
_mkparent(path) = mkpath(dirname(path))


function compute_and_cache_group_metrics_per_df!(
        year::Int,
        profiles_year::OrderedDict{String,<:Any},
        cfg;
        dir::AbstractString = GROUP_DIR,
        overwrite::Bool     = false,
        verbose::Bool       = true)

    for (scen, m_map) in profiles_year
        for (m, slice) in m_map                 # slice :: ProfilesSlice
            variants   = collect(keys(slice.paths))
            n_rep_max  = maximum(length(slice.paths[v]) for v in variants)

            for dem in cfg.demographics
                dem_sym = Symbol(dem)
                pbar = pm.Progress(n_rep_max;
                                desc="[$year|$scen|m=$m|$dem_sym]",
                                barlen=28)

                for rep in 1:n_rep_max
                    cache_path = _perdf(dir, year, scen, m, dem_sym, rep)
                    if isfile(cache_path) && !overwrite
                        pm.next!(pbar); continue
                    end

                    var_map = Dict{Symbol,Vector{DataFrame}}()
                    for var in variants
                        length(slice.paths[var]) < rep && continue
                        df = slice[var, rep]                 # load
                        decode_profile_column!(df)
                        var_map[var] = [df]
                    end
                    isempty(var_map) && (pm.next!(pbar); continue)

                    res = bootstrap_group_metrics(var_map, dem_sym)

                    _mkparent(cache_path)
                    @save cache_path res              # ■ write to disk
                    pm.next!(pbar)
                end
                pm.finish!(pbar); GC.gc()
            end
        end
    end
end

# ────────────────── pass 2: aggregate caches ──────────────────
function accumulate_cached_group_metrics_for_year!(
        year::Int,
        profiles_year::OrderedDict{String,<:Any},
        cfg;
        dir::AbstractString = GROUP_DIR,
        verbose::Bool       = true)

    out = OrderedDict()

    for (scen, m_map) in profiles_year
        scen_out = OrderedDict()
        for (m, slice) in m_map
            variants   = collect(keys(slice.paths))
            n_rep_max  = maximum(length(slice.paths[v]) for v in variants)

            dem_out = OrderedDict()
            for dem in cfg.demographics
                dem_sym = Symbol(dem)
                verbose && @info "  → aggregating $scen, m=$m, dem=$dem_sym"

                accum = Dict{Symbol,Dict{Symbol,Vector{Float64}}}()
                pbar  = pm.Progress(n_rep_max;
                                 desc="[$scen|m=$m|$dem_sym]", barlen=28)

                for rep in 1:n_rep_max
                    cache_path = _perdf(dir, year, scen, m, dem_sym, rep)
                    isfile(cache_path) || error("Missing cache $cache_path")
                    res = nothing; @load cache_path res
                    update_accum!(accum, res, variants)
                    pm.next!(pbar)
                end
                pm.finish!(pbar)
                dem_out[dem_sym] = accum
            end
            scen_out[m] = dem_out
        end
        out[scen] = scen_out
    end
    return out           # scenario ⇒ m ⇒ dem ⇒ metric ⇒ variant ⇒ Vector
end

# ────────────────── public API (drop‑in) ──────────────────
function save_or_load_group_metrics_for_year(year::Int,
                                             profiles_year,
                                             f3_entry;
                                             dir::AbstractString = GROUP_DIR,
                                             overwrite::Bool     = false,
                                             two_pass::Bool      = false,
                                             verbose::Bool       = true)

    final_path = joinpath(dir, "group_metrics_$(year).jld2")

    if isfile(final_path) && !overwrite && !two_pass
        verbose && @warn "Group metrics for $year already cached; loading."
        metrics = nothing; @load final_path metrics; return metrics
    end

    if two_pass
        verbose && @info "Pass 1: computing & caching per‑DataFrame metrics…"
        compute_and_cache_group_metrics_per_df!(year, profiles_year, f3_entry.cfg;
                                                dir=dir, overwrite=overwrite,
                                                verbose=verbose)

        verbose && @info "Pass 2: aggregating cached metrics…"
        metrics = accumulate_cached_group_metrics_for_year!(year, profiles_year,
                                                            f3_entry.cfg;
                                                            dir=dir,
                                                            verbose=verbose)
    else
        verbose && @info "Computing group metrics for year $year (one‑pass)…"
        metrics = apply_group_metrics_for_year_streaming(profiles_year,
                                                         f3_entry.cfg)
    end

    @save final_path metrics
    verbose && @info "Saved aggregated metrics for $year → $final_path"
    return metrics
end



"""
    audit_profiles_year(year; dir = PROFILES_DATA_DIR)

Try reading every encoded profile DF referenced from
`profiles_index_YEAR.jld2`. Returns a Vector of NamedTuples with
(problem=:missing_file|:missing_var|:read_error, and file/scen/m/var/rep/err).
"""
function audit_profiles_year(year::Int; dir::AbstractString = PROFILES_DATA_DIR)
    idxfile = joinpath(dir, "profiles_index_$(year).jld2")
    isfile(idxfile) || error("Index not found: $idxfile")
    profiles_year = JLD2.load(idxfile, "result")  # scenario ⇒ m ⇒ ProfilesSlice

    bad = NamedTuple[]
    # count how many to show a progress bar
    total = sum(length(ps.paths[var]) for (_, m_map) in profiles_year
                                   for (_, ps) in m_map
                                   for var in keys(ps.paths))
    prog = pm.Progress(total; desc = "[audit $year]", barlen = 30)

    for (scen, m_map) in profiles_year
        for (m, ps) in m_map
            for var in keys(ps.paths)
                for (i, fpath) in enumerate(ps.paths[var])
                    if !isfile(fpath)
                        push!(bad, (problem=:missing_file, file=fpath, scenario=scen,
                                    m=m, variant=var, rep=i, err=nothing))
                        pm.next!(prog); continue
                    end
                    # Open without reconstructing the whole DF to check presence of "df"
                    ok, err = true, nothing
                    try
                        JLD2.jldopen(fpath, "r") do f
                            if !haskey(f, "df")
                                ok = false
                                err = "missing variable 'df'"
                            else
                                # force a real read to trigger any reconstruction errors
                                _ = read(f, "df")
                            end
                        end
                    catch e
                        ok = false; err = e
                    end
                    if !ok
                        push!(bad, (problem = err == "missing variable 'df'" ? :missing_var : :read_error,
                                    file=fpath, scenario=scen, m=m, variant=var, rep=i, err=err))
                    end
                    pm.next!(prog)
                end
            end
        end
    end
    pm.finish!(prog)
    return bad
end


"""
    rebuild_profile_file!(ps, iy, cfg, var, rep; dir = PROFILES_DATA_DIR)

Recomputes an encoded profile DF for one (variant, rep) and overwrites it.
Requires: `ps::ProfilesSlice`, `iy::ImputedYear`, `cfg::ElectionConfig`.
"""
function rebuild_profile_file!(ps::ProfilesSlice,
                               iy::ImputedYear,
                               cfg,
                               var::Symbol,
                               rep::Int;
                               dir::AbstractString = PROFILES_DATA_DIR)
    df_imp = iy[var, rep]  # loads imputed replicate
    df = profile_dataframe(df_imp;
                           score_cols = ps.cand_list,  # or ps.cand_syms in your code
                           demo_cols  = cfg.demographics)
    compress_rank_column!(df, ps.cand_list; col = :profile)
    metadata!(df, "candidates", ps.cand_list)
    fprof = ps.paths[var][rep]
    JLD2.@save fprof df
    return fprof
end

"""
    repair_bad_profiles!(year, bad, profiles_year, iy, cfg)

For each entry in `bad` returned by `audit_profiles_year`, rebuild the file.
Skips :missing_file/:missing_var/:read_error uniformly.
"""
function repair_bad_profiles!(year::Int,
                              bad::Vector{<:NamedTuple},
                              profiles_year,
                              iy::ImputedYear,
                              cfg)
    for b in bad
        scen = b.scenario; m = b.m
        ps = profiles_year[scen][m]
        @info "Rebuilding $(basename(b.file))  (scen=$(scen), m=$(m), var=$(b.variant), rep=$(b.rep))"
        rebuild_profile_file!(ps, iy, cfg, b.variant, b.rep)
    end
    return nothing
end



function plot_scenario_year(
    year,
    scenario,
    f3,
    all_meas;
    variant = "mice",
    palette = Makie.wong_colors(),
    figsize = (500,400),
)
    # lookups
    f3_entry   = f3[year]
    cfg        = f3_entry.cfg

    reps_raw   = f3_entry.data
    year_meas  = all_meas[year]
    meas_map   = year_meas[scenario]
    scen_obj   = findfirst(s->s.name==scenario, cfg.scenarios)

    # recompute full candidate set exactly as in generate_profiles_for_year
    sets = unique(map(df->
        compute_candidate_set(df;
            candidate_cols = cfg.candidates,
            m              = cfg.max_candidates,
            force_include  = cfg.scenarios[scen_obj].candidates),
        reps_raw))
    if length(sets)!=1
        msg = "Year $year scenario $scenario: found $(length(sets)) distinct candidate sets; using first."
        @warn msg
        open(CANDLOG,"a") do io
            println(io,"[$(Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))] $msg")
        end
    end
    full_list = sets[1]
    candidate_label = describe_candidate_set(full_list)

    # delegate to your Makie helper
    fig = lines_alt_by_variant(
        meas_map;
        variants        = [variant],
        palette         = palette,
        figsize         = figsize,
        year            = year,
        candidate_label = candidate_label,
    )
    return fig
end



# ──────────────────────────────────────────────────────────────────
# helper: replicate the candidate-set logic from generate_profiles…
# ──────────────────────────────────────────────────────────────────
function _full_candidate_list(cfg, reps_raw, scen_obj)
    sets = unique(map(df ->
        compute_candidate_set(df;
            candidate_cols = cfg.candidates,
            m              = cfg.max_candidates,
            force_include  = scen_obj.candidates),
        reps_raw))
    length(sets) != 1 && @warn "Multiple candidate sets; using first"
    return sets[1]
end



#= # ──────────────────────────────────────────────────────────────────
"""
    plot_group_demographics(
        year::Int,
        scenario::String,
        f3,              # Dict year ⇒ (data,cfg,…)
        all_gm,          # Dict year ⇒ scenario ⇒ m ⇒ dem ⇒ variant ⇒ measure ⇒ Vector
        variant::Symbol  = :mice;
        measures::Vector{Symbol} = [:C, :D, :G],
        figsize::Tuple{Int,Int}  = (1200, 800),
    ) → Figure

Draws a grid (3 × ⌈n/3⌉) of panels — one per demographic — each showing
the three group-metrics C, D, G with 50 % & 95 % ribbons under the chosen
`variant`.  A two-line header gives year, bootstrap count, m-range, and
the human-readable candidate set.  A single legend (rightmost column)
lists only the measures.

Assumes `Makie.wong_colors()` is available (Color 1 → C, 3 → D, 2 → G).
"""
function plot_group_demographics(
    year::Int,
    scenario::String,
    f3,
    all_gm;
    variant ::Symbol       = :mice,
    measures::Vector{Symbol}= [:C, :D, :G],
    figsize ::Tuple{Int,Int}= (1200,800),
)
    # ── look-ups ─────────────────────────────────────────────────────────
    f3ent   = f3[year]
    cfg     = f3ent.cfg
    reps    = f3ent.data
    scenobj = only(filter(s->s.name==scenario, cfg.scenarios))

    gm      = all_gm[year][scenario]                      # m ⇒ dem ⇒ …
    m_vals  = sort(collect(keys(gm)))                     # Vector{Int}
    dems    = cfg.demographics                            # Vector{String}
    nd      = length(dems)

    # ── header info ─────────────────────────────────────────────────────
    cand_lbl = describe_candidate_set(_full_candidate_list(cfg, reps, scenobj))

    
    n_boot     = f3ent.cfg.n_bootstrap

    header1 = "Year = $year • $n_boot bootstraps • m = $(first(m_vals)) … $(last(m_vals))"
    header2 = "$cand_lbl"

    # ── grid layout dims (3×…) ──────────────────────────────────────────
    ncols  = 3
    nrows  = ceil(Int, nd / ncols)
    header_rows = 2          # we occupy first two rows with labels

    # ── colour per MEASURE (same in every panel) ───────────────────────
    wong = Makie.wong_colors()
    palette = Dict(:C => wong[1], :D => wong[3], :G => wong[2])

    # ── create Figure & headers ────────────────────────────────────────
    fig = Figure(resolution = figsize)
    fig[1, 1:ncols] = Label(fig, header1; fontsize=18, halign=:left)
    fig[2, 1:ncols] = Label(fig, header2; fontsize=14, halign=:left)

    # collect first Lines object per measure for the legend
    legend_handles = Dict{Symbol, Lines}()

    xs = Float32.(m_vals)

    for (idx, dem) in enumerate(dems)
        r = header_rows + fld(idx-1, ncols) + 1
        c = mod(idx-1, ncols) + 1

        ax = Axis(fig[r, c];
                  title   = dem,
                  xlabel  = "number of alternatives",
                  ylabel  = idx ≤ ncols ? "value" : "")

        for meas in measures
            # bootstrap arrays along m
            arrays = [ begin
                         vm = gm[m][Symbol(dem)]
                         vals = meas == :G ?
                                (sqrt.(vm[variant][:C] .* vm[variant][:D])) : vm[variant][meas]
                         Float32.(collect(vals))
                       end for m in m_vals ]

            μ   = Float32.(mean.(arrays))
            p25 = Float32.(quantile.(arrays, 0.25f0))
            p75 = Float32.(quantile.(arrays, 0.75f0))
            p05 = Float32.(quantile.(arrays, 0.05f0))
            p95 = Float32.(quantile.(arrays, 0.95f0))

            col = palette[meas]
            band!(ax, xs, p05, p95; color=(col,0.12), linewidth=0)   # 95 %
            band!(ax, xs, p25, p75; color=(col,0.25), linewidth=0)   # 50 %
            ln = lines!(ax, xs, μ; color=col, linestyle=:dot, linewidth=2)

            legend_handles[meas] = get!(legend_handles, meas, ln)
        end
    end

    # ── one legend in the rightmost column ─────────────────────────────
    Legend(fig[header_rows+1:header_rows+nrows, ncols+1],
           [legend_handles[m] for m in measures],
           ["$m • $(variant)" for m in measures];
           orientation = :vertical,
           tellheight  = false)

    return fig
end 
 =#
"""
    plot_group_demographics_lines(
        all_gm, f3, year, scenario;
        variants  = [:zero, :random, :mice],
        measures  = [:C, :D, :G],
        maxcols   = 3,
        n_yticks  = 5,
        palette   = Makie.wong_colors(),
) → Figure

One panel per demographic; **x = number of alternatives (m)**.
For every *(measure, variant)* pair the panel shows

* a translucent band between Q25 and Q75
* a line for the mean.

The title and candidate list match the original
`plot_group_demographics`, but the layout and styling come from
`lines_group_measures_over_m`.
"""
function plot_group_demographics_lines(
        all_gm,
        f3,
        year::Int,
        scenario::String;
        variants      = [:zero, :random, :mice],
        measures      = [:C, :D, :G],
        maxcols::Int  = 3,
        n_yticks::Int = 5,
        palette       = Makie.wong_colors(), clist_size = 60,
)

    # ── data slice & metadata ──────────────────────────────────────────
    gm            = all_gm[year][scenario]                  # m ⇒ (dem ⇒ …)
    m_values_int  = sort(collect(keys(gm)))                 # Vector{Int}
    xs_m          = Float32.(m_values_int)                  # Makie prefers Float32
    demographics  = f3[year].cfg.demographics
    n_demo        = length(demographics)

    scenobj = only(filter(s->s.name==scenario, f3[year].cfg.scenarios))
    cand_lbl = describe_candidate_set(
                 _full_candidate_list(f3[year].cfg, f3[year].data, scenobj))

    # any slice gives the bootstrap length
    sample_slice = gm[first(m_values_int)][Symbol(first(demographics))]
    n_boot       = length(sample_slice[variants[1]][:C])

    # ── colour / style dictionaries ────────────────────────────────────
    measure_cols   = Dict(measures[i] => palette[i] for i in eachindex(measures))
    variant_styles = Dict(:zero => :solid, :random => :dash, :mice => :dot)

    # ── figure geometry ────────────────────────────────────────────────
    ncol       = min(maxcols, n_demo)
    nrow       = ceil(Int, n_demo / ncol)
    title_txt  = "Year $(year) • $(n_boot) bootstraps • m = $(first(m_values_int)) … $(last(m_values_int))"
    fig_width  = max(300*ncol, 10*length(title_txt) + 60)   # widen if title is long
    fig_height = 300*nrow

    fig = Figure(resolution = (fig_width, fig_height))
    rowgap!(fig.layout, 24);  colgap!(fig.layout, 24)

    # headers
    fig[1, 1:ncol] = Label(fig, title_txt;  fontsize = 20, halign = :left)
    fig[2, 1:ncol] = Label(fig,  join(TextWrap.wrap("$cand_lbl"; width=clist_size)); fontsize = 14, halign = :left)
    header_rows = 2

    # legend collectors
    legend_handles = Any[]; legend_labels = String[]

    # ── main panels ────────────────────────────────────────────────────
    for (idx, demo) in enumerate(demographics)
        r, c = fldmod1(idx, ncol)
        ax = Axis(fig[r + header_rows, c];
                  title  = demo,
                  xlabel = "number of alternatives",
                  ylabel = "value",
                  xticks = (xs_m, string.(m_values_int)))            # avoid stretch

        allvals = Float32[]                        # for nice y-ticks

        for meas in measures, var in variants
            vals_per_m = map(m_values_int) do m
                v = gm[m][Symbol(demo)][var]
                arr = meas === :G ? sqrt.(v[:C] .* v[:D]) : v[meas]
                Float32.(arr)
            end

            append!(allvals, vcat(vals_per_m...))

            meds32 = Float32.(mean.(vals_per_m))
            q25s32 = Float32.(map(x -> quantile(x, 0.25f0), vals_per_m))
            q75s32 = Float32.(map(x -> quantile(x, 0.75f0), vals_per_m))

            col = measure_cols[meas]
            sty = variant_styles[var]

            band!(ax, xs_m, q25s32, q75s32; color = (col, 0.20), linewidth = 0)
            ln = lines!(ax, xs_m, meds32;        color = col, linestyle = sty, linewidth = 2)

            if idx == 1
                push!(legend_handles, ln)
                push!(legend_labels, "$(meas) • $(var)")
            end
        end

        # tidy y-ticks
        y_min, y_max = extrema(allvals)
        ticks  = collect(range(y_min, y_max; length = n_yticks))
        ax.yticks[] = (ticks, string.(round.(ticks; digits = 3)))
    end

    # ── legend column ──────────────────────────────────────────────────
    Legend(fig[header_rows+1 : header_rows+nrow, ncol+1],
           legend_handles, legend_labels; tellheight = false)

    # now that the legend is placed, column ncol+1 exists → shrink it
    colsize!(fig.layout, ncol + 1, Relative(0.25))

    resize_to_layout!(fig)
    return fig
end




function compare_demographic_across_scenarios(
        all_gm,
        f3,
        scenario_vec::Vector{Tuple{Int,String}};
        demographic::String,
        variant::Symbol           = :mice,
        measures::Vector{Symbol}  = [:C, :D, :G],
        palette::Vector           = Makie.wong_colors()[1:3],
        n_yticks::Int             = 5,
        base_width::Int           = 400,
        base_height::Int          = 360
)

    # ── colour helper (lighten / darken without unqualified imports) ──
    ΔL = Dict(:C=>0, :D=>+20, :G=>-15)
    function shade(rgb::Colors.RGB, meas)
        lch = convert(Colors.LCHab, rgb)
        newL = clamp(lch.l + ΔL[meas], 0, 100)
        convert(Colors.RGB, Colors.LCHab(newL, lch.c, lch.h))
    end

    base_rgbs = Colors.RGB.(palette[1:length(measures)])
    measure_cols = Dict(measures[i] => shade(base_rgbs[i], measures[i])
                        for i in eachindex(measures))

    # ── slice metadata from the first scenario ────────────────────────
    y0, s0   = scenario_vec[1]
    gm0      = all_gm[y0][s0]
    m_vals   = sort(collect(keys(gm0)))           # Int
    xs_m     = Float32.(m_vals)
    n_panels = length(scenario_vec)
    demo_sym = Symbol(demographic)

    n_boot = length(gm0[first(m_vals)][Symbol(demographic)][variant][:C])

    # helper for candidate label
    candidate_label(y, s) = begin
        cfg  = f3[y].cfg
        scenobj = only(filter(t -> t.name == s, cfg.scenarios))
        describe_candidate_set(_full_candidate_list(cfg, f3[y].data, scenobj))
    end

    # ── global y-limits ───────────────────────────────────────────────
    global_all = Float32[]
    for (yr, sc) in scenario_vec, meas in measures, m in m_vals
        v = all_gm[yr][sc][m][Symbol(demographic)][variant]
        arr = meas === :G ? sqrt.(v[:C] .* v[:D]) : v[meas]
        append!(global_all, Float32.(arr))
    end
    y_min, y_max = extrema(global_all)
    y_ticks = collect(range(y_min, y_max; length = n_yticks))

    # ── figure scaffold (3 rows) ──────────────────────────────────────
    fig = Figure(resolution = (base_width*n_panels, base_height),
                 layout = (3, n_panels))
    rowgap!(fig.layout, 20); colgap!(fig.layout, 30)

    header_txt = "$demographic • number of alternatives = $(first(m_vals))…$(last(m_vals)) • $n_boot bootstraps"
    fig[1, 1:n_panels] = Label(fig, header_txt; fontsize = 22, halign = :center)

    legend_handles = Lines[]; legend_labels = String[]

    # ── panel loop ────────────────────────────────────────────────────
    for (i, (yr, sc)) in enumerate(scenario_vec)
        gm_slice = all_gm[yr][sc]
        wrapped_title = join(TextWrap.wrap("Year $yr — $(candidate_label(yr, sc))"; width = 50))

        ax = Axis(fig[2, i];
                  title     = wrapped_title,
                    titlefont = "sans",
                  titlesize = 14,
                  titlegap  = 8,
                  xlabel    = "number of alternatives",
                  ylabel    = "value",
                  xticks    = (xs_m, string.(m_vals)),
                  limits    = (nothing, (y_min, y_max)),
                  yticks    = (y_ticks, string.(round.(y_ticks; digits=2))))

        for meas in measures
            col = measure_cols[meas]
            meds = Float32[]; q25s = Float32[]; q75s = Float32[]; p05s = Float32[]; p95s = Float32[]

            for m in m_vals
                v = gm_slice[m][Symbol(demographic)][variant]
                vals32 = Float32.(meas === :G ? sqrt.(v[:C] .* v[:D]) : v[meas])

                push!(meds, median(vals32))
                push!(q25s, quantile(vals32, 0.25f0)); push!(q75s, quantile(vals32, 0.75f0))
                push!(p05s, quantile(vals32, 0.05f0)); push!(p95s, quantile(vals32, 0.95f0))
            end

            band!(ax, xs_m, p05s, p95s; color = (col, 0.12), linewidth = 0)
            band!(ax, xs_m, q25s, q75s; color = (col, 0.25), linewidth = 0)
            ln = lines!(ax, xs_m, meds; color = col, linewidth = 2)

            if i == 1
                push!(legend_handles, ln)
                push!(legend_labels, string(meas))
            end
        end
    end

    fig[3, 1:n_panels] = Legend(fig, legend_handles, legend_labels;
                                orientation = :horizontal,
                                framevisible = false,
                                halign = :center)

    resize_to_layout!(fig)
    return fig
end





"""
    save_plot(fig, year, scenario, cfg; variant, dir = "imgs", ext = ".png")

Save `fig` under `dir/`, creating the directory if needed.  
The file name pattern is:

    {year}_{scenario}_{variant}_B{n_bootstrap}_M{max_m}_{yyyymmdd-HHMMSS}{ext}
"""
function save_plot(fig, year::Int, scenario::AbstractString, cfg;
                   variant::AbstractString,
                   dir::AbstractString = "imgs",
                   ext::AbstractString = ".png")

    # 1. make sure the directory exists
    mkpath(dir)

    # 2. assemble a human-readable, reproducible file name
    time_stamp = Dates.format(now(), "yyyymmdd-HHMMSS")
    max_m      = maximum(cfg.m_values_range)
    fname      = joinpath(dir,
        string(year, '_', scenario, '_', variant,
               "_B", cfg.n_bootstrap,
               "_M", max_m,
               '_', time_stamp, ext))

    # 3. save
    save(fname, fig; px_per_unit = 2)
    @info "saved plot → $fname"
    return fname
end

"""
    median_neff(all_meas, year, scenario;
                variant    = :mice,
                m_max      = typemax(Int),
                aggregate  = false)

Compute the median effective number of reversal pairs  
N_eff = 1 / calc_reversal_HHI   for the requested
`year` / `scenario` / `variant`.

* If `aggregate == false` (default) → return `Dict(m ⇒ median N_eff(m))`.
* If `aggregate == true`            → return a single `Float64` with the
  median taken over *all m ≤ m_max*.

`variant` may be a `Symbol` or `String`.  Set `m_max` to 5 if you want to
exclude larger numbers of alternatives.
"""
function median_neff(all_meas,
                     year::Integer,
                     scenario::AbstractString;
                     variant    = :mice,
                     m_max::Int = typemax(Int),
                     aggregate::Bool = false)

    var = Symbol(variant)                       # normalise the key

    scen_map = all_meas[year][scenario]         # m ⇒ measure ⇒ …

    # collect medians per m
    med_per_m = Dict{Int,Float64}()
    for (m, mdict) in scen_map
        m > m_max && continue
        hhi_vec = mdict[:calc_reversal_HHI][var]        # already normalised
        med_per_m[m] = median(1.0 ./ hhi_vec)           # N_eff = 1/HHI_*
    end

    return aggregate ? median(values(med_per_m)) : med_per_m
end



function enrp_table(all_meas, f3,
year::Integer, scenario::AbstractString;
variant = :mice, m_max::Int = typemax(Int))


var   = Symbol(variant)
data  = all_meas[year][scenario]              # m ⇒ measure ⇒ …
m_vec = sort([m for m in keys(data) if m ≤ m_max])

tbl = OrderedCollections.OrderedDict{Int,NamedTuple}()

for m in m_vec
    hhi_vec = data[m][:calc_reversal_HHI][var]    # already normalised
    enrp    = median(1.0 ./ hhi_vec)
    tbl[m]  = (enrp = enrp, max = factorial(big(m)) / 2)  # big() avoids overflow
end
return tbl
end





"""
    hhi_table(all_meas, f3, year, scenario;
              variant = :mice,
              m_max   = typemax(Int))

Return an `OrderedDict m ⇒ (hhi = median HHI_*, min = 2/factorial(m))`
sorted by `m`.  Works just like `enrp_table`.
"""
function hhi_table(all_meas,
                   f3,
                   year::Integer,
                   scenario::AbstractString;
                   variant = :mice,
                   m_max::Int = typemax(Int))

    var   = Symbol(variant)
    data  = all_meas[year][scenario]                 # m ⇒ measure ⇒ …
    m_vec = sort([m for m in keys(data) if m ≤ m_max])

    tbl = OrderedCollections.OrderedDict{Int,NamedTuple}()

    for m in m_vec
        hhi_vec = data[m][:calc_reversal_HHI][var]   # already normalised
        tbl[m]  = (
            hhi = median(hhi_vec),
            min = 2.0 / float(factorial(big(m))),    # 1 / (#pairs)
        )
    end
    return tbl
end



function polar_table(all_meas,
                     _f3,                         # unused, kept for API parity
                     year::Integer,
                     scenario::AbstractString;
                     variant    = :mice,
                     m_max::Int = typemax(Int))

    var   = Symbol(variant)
    data  = all_meas[year][scenario]                   # m ⇒ measure ⇒ …
    m_vec = sort([m for m in keys(data) if m ≤ m_max])

    tbl = OrderedDict{Int,NamedTuple}()

    for m in m_vec
        md  = data[m]

        psi_vec  = md[Symbol("Ψ")][var]
        rhhi_vec = md[:fast_reversal_geometric][var]
        r_vec    = md[:calc_total_reversal_component][var]

        tbl[m] = (
            psi  = median(psi_vec),
            rhhi = median(rhhi_vec),
            r    = median(r_vec),
        )
    end
    return tbl
end
