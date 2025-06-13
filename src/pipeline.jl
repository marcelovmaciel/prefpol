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


