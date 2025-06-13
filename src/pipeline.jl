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
