
# === Project and Data Setup ===

function load_project_path()
    return dirname(Pkg.project().path)
end

function load_spss_file(path::String)
    RCall.reval("library(haven)")
    RCall.reval("library(mice)")

    @rput path  # this macro works fine for variables
    return rcopy(R"read_sav(path)")
end



# === Score Processing ===
function build_candidate_score_distributions(df::DataFrame, candidatos::Vector{String})
    return Dict(c => countmap(df[!, c]) for c in candidatos)
end

#= function extract_unique_scores(countmaps::Dict{String,<:AbstractDict})
   
    return sort!(unique(vcat([keys(cm) for cm in values(countmaps)]...)))
end
 =#

function extract_unique_scores(countmaps::Dict{String,<:AbstractDict})
    all_keys = Int[]
    for cm in values(countmaps)
        append!(all_keys, Int.(keys(cm)))
    end
    return sort!(unique(all_keys))
end


# function convert_keys_to_int(dict)
#     return Dict(Int(k) => v for (k, v) in dict)
# end

function convert_keys_to_int(dict)
    return Dict(
        if k isa Integer
            k
        elseif k isa AbstractFloat
            r = round(Int, k)
            if !isapprox(k, r; atol=1e-8)
                throw(ArgumentError("Expected near-integer key, got $k"))
            end
            r
        else
            throw(ArgumentError("Unsupported key type: $(typeof(k))"))
        end => v
        for (k, v) in dict
    )
end


function sanitize_countmaps(countmaps::Dict{String,<:AbstractDict})
    return Dict(c => convert_keys_to_int(cm) for (c, cm) in countmaps)
end


function compute_dont_know_her(countmaps::Dict{String,Dict{Int,Int}}, nrespondents::Int)
    return sort([
        (k, 100 * sum(get(v, code, 0) for code in (96, 97, 98,99)) / nrespondents)
        for (k, v) in countmaps
    ], by = x -> x[2])
end



function prepare_scores_for_imputation_int(df::DataFrame,
    score_cols::Vector{String};
    extra_cols::Vector{String}=String[])
    # (1) Split truly numeric from anything that isn't
    numeric_cols   = Base.filter(c -> eltype(df[!, c]) <: Union{Missing, Real}, score_cols)
    nonnumeric     = setdiff(score_cols, numeric_cols)
    if !isempty(nonnumeric)
        @warn "prepare_scores_for_imputation_int: skipping non‑numeric columns $(nonnumeric)"
    end

    # (2) Work only on the numeric score columns
    scores_int = mapcols(col -> Int.(col), df[:, numeric_cols])
    declared   = Impute.declaremissings(scores_int; values = (96, 97, 98, 99))

    # (3) Append any extra (demographic) columns, untouched
    return isempty(extra_cols) ? declared : hcat(declared, df[:, extra_cols])
end

#= function prepare_scores_for_imputation_categorical(df::DataFrame, cols::Vector{String})
    declared = prepare_scores_for_imputation_int(df, cols)
    return mapcols(x -> categorical(x, ordered=true), declared)
end
 =#

"""
    prepare_scores_for_imputation_categorical(df, score_cols; extra_cols = String[])

Same idea as the `_int` version but returns the scores as **ordered categoricals**.
Demographics are still appended unchanged.
"""
function prepare_scores_for_imputation_categorical(df::DataFrame,
                                                   score_cols::Vector{String};
                                                   extra_cols::Vector{String}=String[])
    declared = prepare_scores_for_imputation_int(df, score_cols; extra_cols = String[])
    declared_cat = mapcols(col -> categorical(col, ordered = true), declared)
    return isempty(extra_cols) ? declared_cat : hcat(declared_cat, df[:, extra_cols])
end


# === Slice Top  ===
function get_most_known_candidates(dont_know_her::Vector{Tuple{String, Float64}}, how_many)
    most_known_candidates = [x[1] for x in dont_know_her[1:how_many]]    
    return most_known_candidates
end

#= """
    select_top_candidates(countmaps::Dict, nrespondents::Int;
                          m::Int,
                          force_include::Vector{String}=String[])

Returns a vector of length **≤ m** with

1. every name in `force_include` (duplicates ignored ➔ unique list)
2. the remaining slots filled by the lowest “don't-know-her” share.

Sanity checks:
* if `length(force_include) > m`  → `ArgumentError`
* if after filling we still have < m distinct names (e.g. survey listed
  only k < m candidates) we emit a `@warn`.

The returned vector preserves the order:
`force_include` first (as given) then the popularity-sorted remainder.
"""
function select_top_candidates(countmaps::Dict{String,<:AbstractDict},
                               nrespondents::Int;
                               m::Int,
                               force_include::Vector{String}=String[])

    # ——— sanity -------------------------------------------------------------
    inc = unique(force_include)
    if length(inc) > m
        throw(ArgumentError("""
        force_include has $(length(inc)) names but m = $m;
        drop some names or request a larger m."""))
    end

    # ——— 1. names that must appear -----------------------------------------
    selected = copy(inc)                         # preserve caller order
    set_inc  = Set(inc)                          # for fast lookup

    # ——— 2. the usual popularity ranking (excluding forced ones) ------------
    dont_know = compute_dont_know_her(countmaps, nrespondents)  # [(name,pct), …]
    for (name, _) in dont_know
        name ∈ set_inc && continue               # already selected
        push!(selected, name)
        length(selected) == m && break
    end

    if length(selected) < m
        @warn "Only $(length(selected)) unique candidates found; requested $m."
    end
    return selected
end =#

#= 
function select_top_candidates(countmaps::Dict{String,<:AbstractDict},
                               nrespondents::Int;
                               m::Int,
                               force_include::Vector{String}=String[])

    inc = unique(force_include)
    if length(inc) > m
        throw(ArgumentError("force_include has $(length(inc)) names but m = $m"))
    end

    # ---------- popularity list, already sorted ascending by “don't-know-her” ----
    poplist = [name
               for (name, _) in compute_dont_know_her(countmaps, nrespondents)
               if name ∉ inc]                         # drop forced names, duplicates

    needed  = m - length(inc)
    extra   = poplist[1:needed]                       # may be empty

    selected = vcat(inc, extra)

    if length(selected) < m
        @warn "Only $(length(selected)) unique candidates available; requested $m."
    end
    return selected
end =#

function select_top_candidates(countmaps::Dict{String,<:AbstractDict},
                               nrespondents::Int;
                               m::Int,
                               force_include::Vector{String}=String[])

    # remove duplicates, keep order
    inc = unique(force_include)

    # truncate if too many
    if length(inc) > m
        @warn "force_include has more than $m names; truncating to first $m."
        inc = inc[1:m]
    end

    # ---------- popularity list, already sorted ascending by “don't-know-her” ----
    poplist = [name
               for (name, _) in compute_dont_know_her(countmaps, nrespondents)
               if name ∉ inc]                         # drop forced names, duplicates

    needed  = m - length(inc)
    extra   = needed > 0 ? poplist[1:min(needed, length(poplist))] : String[]

    selected = vcat(inc, extra)

    if length(selected) < m
        @warn "Only $(length(selected)) unique candidates available; requested $m."
    end

    return selected
end


function compute_candidate_set(scores_df::DataFrame;
                               candidate_cols,
                               m::Int,
                               force_include::Vector{String} = String[])

    countmaps     = build_candidate_score_distributions(scores_df, candidate_cols)
    countmaps2    = sanitize_countmaps(countmaps)
    nrespondents  = nrow(scores_df)

    return select_top_candidates(countmaps2, nrespondents;
                                 m = m, force_include = force_include)
end



# filter scores_df to just have the most known candidates

function get_df_just_top_candidates(df::DataFrame, how_top::Int; demographics = String[] )
    most_known_candidates = get_most_known_candidates(dont_know_her, how_top)
    return df[!, vcat(most_known_candidates, demographics)]
end

function get_df_just_top_candidates(df::DataFrame, which_ones; demographics = String[]) 
    return df[!, vcat(which_ones, demographics)]
end



#= 

# note the m=1 here! I could play with stats of the imputation....
# another approach would be to boostrap the original dataset THEN apply to EACH of these! 
const GLOBAL_R_IMPUTATION = let
    function f(df::DataFrame)
        seed = rand(1:10^6)
        RCall.reval("set.seed($seed)") 
        R"""
        
        library(mice)
        options(nnet.MaxNWts = 10000L)  
        df <- as.data.frame($df)
        init <- mice(df, maxit = 0, print = FALSE)
        meth <- init$method
        pred <- make.predictorMatrix(df)
        diag(pred) <- 0

        # Safe method assignment
        for (v in names(df)) {
          u <- unique(na.omit(df[[v]]))
          n_cat <- length(u)

          if (n_cat <= 1) {
            meth[v] <- ""                    # Skip constant or all-missing vars
          } else if (is.factor(df[[v]])) {
            # Polyreg is fragile: fallback if >100 predictors or lots of categories
            if (sum(pred[v, ]) > 100 || n_cat > 10) {
              meth[v] <- "sample"
            } else {
              meth[v] <- "polyreg"
            }
          } else if (is.numeric(df[[v]])) {
            meth[v] <- "pmm"
          }
        }

        imp <- mice(df,
                    m = 1,
                    method = meth,
                    predictorMatrix = pred,
                    seed = 123,
                    printFlag = FALSE)

        completed_df <- complete(imp, 1)
        """
        return rcopy(DataFrame, R"completed_df")
    end
end =#


const GLOBAL_R_IMPUTATION = let
    function f(df::DataFrame; m::Int = 1)
        # random seed to keep bootstrap independence
        seed = rand(1:10^6)
        RCall.reval("set.seed($seed)")

        R"""
        suppressPackageStartupMessages(library(mice))

        df <- as.data.frame($df)

        # ---------- boilerplate ----------
        init <- mice(df, maxit = 0, print = FALSE)
        meth <- init$method                   # default methods
        pred <- make.predictorMatrix(df)
        diag(pred) <- 0                       # no self-prediction

        # ---------- customise methods ----------
        for (v in names(df)) {
          col <- df[[v]]
          if (all(is.na(col)) || length(unique(na.omit(col))) <= 1) {
            meth[v] <- ""                     # constant or all-missing
          } else if (is.factor(col)) {
            n_cat <- nlevels(col)
            if (n_cat == 2) {
              meth[v] <- "logreg"             # binomial GLM
            } else {                          # 3+ categories
              meth[v] <- "cart"               # safe, no weight explosion
            }
          } else if (is.numeric(col)) {
            meth[v] <- "pmm"
          }
        }

        # ---------- one imputation ----------
        imp <- mice(df,
                    m               = $m,
                    method          = meth,
                    predictorMatrix = pred,
                    printFlag       = FALSE)

        completed_df <- complete(imp, 1)
        """
        return rcopy(DataFrame, R"completed_df")
    end
end

# === End of Slice Top ===
# === Imputation ===


function imputation_variants(df::DataFrame,
    candidates::Vector{String},
    demographics::Vector{String};
    most_known_candidates::Vector{String}=String[])

# 1 ─ Determine which score columns to use
use_cols = isempty(most_known_candidates) ? candidates : most_known_candidates

# 2 ─ Subset to relevant columns (if top-candidates requested)
df_subset = isempty(most_known_candidates) ? df : get_df_just_top_candidates(df, use_cols; demographics = demographics)

# 3 ─ Prepare imputation tables (same for both branches now)
scores_int  = prepare_scores_for_imputation_int(df_subset, use_cols; extra_cols = demographics)
scores_cat  = prepare_scores_for_imputation_categorical(df_subset, use_cols; extra_cols = demographics)

# 4 ─ Apply imputation variants
imputed0    = Impute.replace(scores_int, values = 0)
imputedRnd  = Impute.impute(scores_cat, Impute.SRS(; rng = MersenneTwister()))
imputedM    = GLOBAL_R_IMPUTATION(scores_cat)

return (zero = imputed0,
random = imputedRnd,
mice = imputedM)
end




# === End of Imputation ===

#= TODO: here is where I should be writing the bootstrap.
I should be using: 
  - the imputed data from mice
  - the imputed data from SRS
  - the 0 score imputation 
  - the original data with missing values  (for the measures in which this is possible!)

Now, I must to remember to use the survey weights in the bootstrap.  
=#

    function weighted_bootstrap(data::DataFrame, weights::Vector{Float64}, B::Int)
    n = nrow(data)
    boot_samples = Vector{DataFrame}(undef, B)
    
    for b in 1:B
        idxs = sample(1:n, Weights(weights), n; replace=true)
        boot_samples[b] = data[idxs, :]
    end
    
    return boot_samples
end

"""
    bootstrap_variants(imps; weights, B = 10, keys = (:zero, :random, :mice))

• `imps`    – NamedTuple returned by `imputed_variants`  
• `weights` – vector of survey weights (same length as nrow of the imputations)  
• `B`       – number of bootstrap replications (default 10)  
• `keys`    – which variants inside `imps` to bootstrap (default skips the
              raw tables with missings)

Returns `Dict{Symbol, Vector{DataFrame}}`.
"""
function bootstrap_variants(imps::NamedTuple;
                            weights::AbstractVector,
                            B::Int = 10,
                            keys::Tuple = (:zero, :random, :mice))

    Dict(k => weighted_bootstrap(imps[k], weights, B) for k in keys)
end


function checked_weighted_bootstrap(df::DataFrame,
                                    code_to_r,
                                    weights::Vector{Float64},
                                    B::Int)

    n = nrow(df)
    boot = Vector{DataFrame}(undef, B)

    # max valid code (for cheap comparison)
    maxcode = maximum(keys(code_to_r))

    for b in 1:B
        idxs = sample(1:n, Weights(weights), n; replace = true)
        sub  = df[idxs, :]

        # inexpensive check: any(code > max?)  — if yes, run full scan
        if any(code -> code > maxcode || !haskey(code_to_r, code),
               sub.profile)
            bad = unique(filter(code -> !haskey(code_to_r, code),
                                sub.profile))
            error("Bootstrap replicate $b contained unknown codes: $bad")
        end
        boot[b] = sub
    end
    return boot
end

function bootstrap_encoded_variants(enc_imps;
weights::AbstractVector,
B::Int = 10,
keys::Tuple = (:zero, :random, :mice))
result = Dict{Symbol,Dict}()

for k in keys
    enc  = enc_imps[k]                          # Dict with :df, :r_to_code, :code_to_r
    df   = enc[:df]
    code = enc[:code_to_r]
    
    boot = checked_weighted_bootstrap(df, code, weights, B)   # Vector{DataFrame}

    result[k] = Dict(
        :dfs       => boot,
        :r_to_code => enc[:r_to_code],
        :code_to_r => enc[:code_to_r],
    )
end
return result
end


function get_row_candidate_score_pairs(row, score_cols) 
    Dict(Symbol(c) => row[c] for c in score_cols)
end

function  get_order_dict(score_dict) 
    unique_scores = sort(unique(values(score_dict)); rev = true)
    lookup = Dict(s => r for (r,s) in enumerate(unique_scores))
    Dict(k => lookup[v] for (k,v) in score_dict)
end

function force_scores_become_linear_rankings(score_dict; rng=MersenneTwister())

    grouped = Dict(score => Symbol[] for score in unique(values(score_dict)))
    
    for (cand, score) in score_dict
        push!(grouped[score], cand)
    end

    sorted_scores = sort(collect(keys(grouped)), rev=true)
    linear_ranking = Dict{Symbol, Int}()
    next_rank = 1

    for score in sorted_scores
        cands = grouped[score]
        shuffle!(rng, cands)
        for cand in cands
            linear_ranking[cand] = next_rank
            next_rank += 1
        end
    end

    return linear_ranking
end




# ————————————————————— profile builders ——————————————————————

function build_profile(df::DataFrame;
                       score_cols::Vector,
                       rng  = Random.GLOBAL_RNG,
                       kind::Symbol = :linear)   # :linear or :weak
    f = kind === :linear ? force_scores_become_linear_rankings : get_order_dict
    score_dicts = map(row -> get_row_candidate_score_pairs(row, score_cols),
                      eachrow(df))
    
    return map(sd -> f(sd), score_dicts)
end


function profile_dataframe(df::DataFrame;
                           score_cols::Vector,
                           demo_cols::Vector,
                           rng  = Random.GLOBAL_RNG,
                           kind::Symbol = :linear)
    prof = build_profile(df; score_cols = score_cols, rng = rng, kind = kind)
    return DataFrame(profile = prof) |> (d -> hcat(d, df[:, demo_cols]))
end




function encode_imputation_variants(variants, score_cols, demographics)
    # pairs(variants) yields (key, value) for NamedTuple, Dict, etc.
    return Dict(k => encoded_imputed_df(v, score_cols, demographics)
                for (k, v) in pairs(variants))
end 





function normalize_scores!(df::DataFrame, score_syms::Vector{Symbol})
    for s in score_syms
        col = df[!, s]

        if eltype(col) <: Union{Missing, Int}
            continue                               # already nice
        end

        allowmissing!(df, s)                       # ensure we can keep missings

        if eltype(col) <: CategoricalValue
            # turn categorical → underlying string/int
            raw = levels(col)[levelcode.(col)]
        else
            raw = collect(col)                     # plain Vector
        end

        df[!, s] = Union{Missing, Int}[           # overwrite with Int / missing
            x === missing           ? missing :
            x isa Int               ? x :
            x isa AbstractString    ? tryparse(Int, x) :
            x isa Real              ? Int(x) :
                                      missing
            for x in raw
        ]
    end
    return df
end 


function normalize_scores_int!(df::DataFrame, score_syms::Vector{Symbol})
    for s in score_syms
        col = df[!, s]

        eltype(col) <: Int && continue  # already good

        any(ismissing, col) &&
            error("Column $s still contains missing after imputation")

        vec = Vector{Int}(undef, length(col))

        @inbounds for i in eachindex(col)
            x = col[i]
            if x isa CategoricalValue
                raw = levels(x)[levelcode(x)]     # stored value (Int or String)
                vec[i] = raw isa Int ? raw : parse(Int, raw)
            else
                vec[i] = x isa Int ? x :
                          x isa AbstractString ? parse(Int, x) :
                          Int(x)                             # Float64 etc.
            end
        end
        df[!, s] = vec
    end
    return df
end

function encoded_imputed_df(df, score_cols, demographics)
    foo = profile_dataframe(normalize_scores!(copy(df), Symbol.(score_cols)); score_cols = score_cols, demo_cols = demographics)
    # foo = profile_dataframe(normalize_scores_int!(copy(df), Symbol.(score_cols)); score_cols = score_cols, demo_cols = demographics)
    # second variant might be faster 
    # Build full codebook on *all* rows
    ranks = unique(foo.profile)                     # every ranking present
    r_to_code = Dict(r => i for (i,r) in enumerate(ranks))
    code_to_r = Dict(i => r for (i,r) in enumerate(ranks))

    # --- encode in-place, checking every row -------------------------------
    foo.profile = [ get(r_to_code, r) do
                        error("Row with unseen ranking slipped through the map.")
                    end for r in foo.profile ]

    return Dict(
        :df        => foo,          # fully encoded DataFrame
        :r_to_code => r_to_code,
        :code_to_r => code_to_r,
    )
end



function decode_profile_vec(prof, code_to_r)
    return map(code -> code_to_r[code], prof)
end




function check_keys(bt, key)
    println("values in profile: ", sort(unique(vcat(map(x-> unique(x.profile),bt[key][:dfs])...))) )
    println("type of value in profile: ", typeof(bt[key][:dfs][1].profile[1]))
    println("original values: " ,sort(Int.(keys(bt[key][:code_to_r]))))
    println("type of value in original: ", typeof(bt[key][:code_to_r]))
end




@inline function dict2svec(d::Dict{Symbol,<:Integer}; cs::Vector{Symbol}=cands,
                           higher_is_better::Bool=false)
    # 1. pack the m scores into an isbits StaticVector
    m = length(cs)                                # number of candidates
    vals = SVector{m,Int}(map(c -> d[c], cs))

    # 2. permutation that sorts those scores
    perm = sortperm(vals; rev = higher_is_better)           # Vector{Int}

    # 3. return as SVector{m,UInt8} (10 B if m ≤ 10)
    return SVector{m,UInt8}(perm)
end




# keep original    pool[code]  behaviour
decode_rank(code::Integer,      pool) = pool[code]

# call is a no-op if you pass the SVector itself
decode_rank(r::SVector, _) = r



function compress_rank_column!(df::DataFrame, cands; col::Symbol=:profile)
    # 1. Dict → SVector
    
    sv = [dict2svec(r[col],cs = cands) for r in eachrow(df)]  # one tiny allocation per row

    # 2. pool identical SVectors (UInt16 index)
    pooled = PooledArray(sv; compress = true)

    # 3. overwrite in-place; let Dict objects be GC’d
    df[!, col] = pooled
    GC.gc()                       # reclaim Dict storage promptly
    return pooled.pool            # decoder lookup table
end



@inline function perm2dict(perm::AbstractVector{<:Integer},
                           cs::Vector{Symbol})
    d = Dict{Symbol,Int}()
    @inbounds for (place, idx) in pairs(perm)          # place = 1,2,…
        d[cs[idx]] = place
    end
    return d
end


perm_to_dict = @inline perm2dict


function decode_profile_column!(df::DataFrame)
    eltype(df.profile) <: Dict && return df            # nothing to do

    cand_syms = metadata(df, "candidates")
    col       = df.profile
    decoded   = Vector{Dict{Symbol,Int}}(undef, length(col))

    if col isa PooledArray
        pool = col.pool
        for j in eachindex(col)
            perm = decode_rank(col[j], pool)
            decoded[j] = perm_to_dict(perm, cand_syms)
        end
    else                                               # plain Vector{SVector}
        for j in eachindex(col)
            decoded[j] = perm_to_dict(col[j], cand_syms)
        end
    end

    df[!, :profile] = decoded
    return df
end


@inline function decode_each!(var_map)
    for vec in values(var_map)          # vec::Vector{DataFrame}
        decode_profile_column!(vec[1])  # length == 1 in streaming path
    end
end
