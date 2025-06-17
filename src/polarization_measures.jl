
function find_reversal_pairs(unique_rankings::Vector{<:NTuple{N, String}}) where {N}
    paired_accum = Vector{Tuple{NTuple{N,String}, Int, NTuple{N,String}, Int}}()
    unpaired_accum = Vector{Tuple{NTuple{N,String}, Int}}()
    paired_indices = Set{Int}()

    for (i, ranking) in enumerate(unique_rankings)
        # Skip if already paired
        if i in paired_indices
            continue
        end

        # Reverse the ranking
        rev_ranking = reverse(ranking)
        found_index = nothing

        # Search for the reversed ranking in the remaining items
        for j in (i+1):length(unique_rankings)
            if j in paired_indices
                continue
            end
            if unique_rankings[j] == rev_ranking
                found_index = j
                break
            end
        end

        if isnothing(found_index)
            # No reversal found – record as unpaired
            push!(unpaired_accum, (ranking, i))
        else
            # Found a reversal – record the pair and mark both indices as paired
            push!(paired_accum, (ranking, i, rev_ranking, found_index))
            push!(paired_indices, i)
            push!(paired_indices, found_index)
        end
    end

    return paired_accum, unpaired_accum
end

function print_reversal_results(paired_accum, unpaired_accum)
    println("Paired Reversals:")
    if isempty(paired_accum)
        println("No reversal pairs found.")
    else
        for (r, i, rev_r, j) in paired_accum
            println("Index $i: $r  <==> Index $j: $rev_r")
        end
    end

    println("\nUnpaired Rankings:")
    if isempty(unpaired_accum)
        println("No unpaired rankings.")
    else
        for (r, i) in unpaired_accum
            println("Index $i: $r")
        end
    end
end

# paired, unpaired = find_reversal_pairs(unique_top_rankings)




"""Return an iterator of local reversal values (2 * min(prop_i, prop_j)) 
   for each reversal pair in `paired_accum`."""
function local_reversal_values(
    paired_accum::Vector{<:Tuple{NTuple{N,String},Int,NTuple{N,String},Int}},
    proportion_rankings::Dict{NTuple{N,String},Float64}
) where N
    return (
        2 * min(proportion_rankings[p[1]], proportion_rankings[p[3]])
        for p in paired_accum
    )
end

"""Sum of local reversal components: Σ (2 * min(prop_i, prop_j))"""
function calc_total_reversal_component(paired_accum, proportion_rankings::Dict)
    sum(local_reversal_values(paired_accum, proportion_rankings))
end



"""Sum of squares of local reversal components: Σ (2 * min(prop_i, prop_j))^2"""
function calc_reversal_HHI(paired_accum, proportion_rankings::Dict)
    sum(x^2 for x in local_reversal_values(paired_accum, proportion_rankings))
end

"""Geometric measure of reversal component: 
   sqrt( (sum of local reversals) * (sum of squares of local reversals) )"""
function calc_reversal_geometric(paired_accum, proportion_rankings::Dict)
    reversal_component = calc_total_reversal_component(paired_accum, proportion_rankings)
    reversal_hhi       = calc_reversal_HHI(paired_accum, proportion_rankings)
    sqrt(reversal_component * reversal_hhi)
end



"""Faster one-pass geometric measure: 
   accumulate sum and sum of squares in one loop, then sqrt(...)"""
function fast_reversal_geometric(paired_accum, proportion_rankings::Dict)
    total_reversal = 0.0
    total_hhi = 0.0
    for x in local_reversal_values(paired_accum, proportion_rankings)
        total_reversal += x
        total_hhi      += x^2
    end
    sqrt(total_reversal * total_hhi)
end



# Count how many rankings rank candidate1 higher than candidate2.
function nab(candidate1, candidate2, profile::Vector{<:Dict})
    return count(ranking -> ranking[candidate1] < ranking[candidate2], profile)
end

# Compute the absolute difference between the counts for candidate1 and candidate2.
function dab(candidate1, candidate2, profile::Vector{<:Dict})
    return abs(nab(candidate1, candidate2, profile) - nab(candidate2, candidate1, profile))
end



function Ψ(p)
    cs = collect(keys(p[1]))
    candidate_pairs = collect(combinations(cs,2))

    m_choose_2 = length(candidate_pairs)
    n = length(p)

    ∑ = sum 

    can_polarization = ∑([(n-dab(c1,c2,p)) for (c1,c2) in candidate_pairs])/(n*m_choose_2)
    return(can_polarization)

end # TODO: note this function can be applied even to incomplete preferences, since only operates on pairs 





function w(c1, c2, ranking)
    1/(ranking[c1] + ranking[c2])
end




"""
    ptilde(c1, c2, profile)

Compute the weighted proportion of voters who prefer `c1` to `c2`:
  
  \tilde{p}_{c1,c2} = (∑ over i of w(c1, c2, i) * 1{c1 ≻_i c2}) / (∑ over i of w(c1, c2, i))

where `w` is a weighting function, and `1{c1 ≻_i c2}` is an indicator
that voter `i` ranks `c1` strictly above `c2`.
"""
function ptilde(c1, c2, profile::Vector{<:Dict{Symbol, Int}})
    numerator   = 0.0
    denominator = 0.0
    for ranking in profile
        local_weight = w(c1, c2, ranking)
        denominator += local_weight
        if ranking[c1] < ranking[c2]
            numerator += local_weight
        end
    end
    return denominator == 0 ? 0.0 : (numerator / denominator)
end


# The weighted Psi measure using 1 - |2 p_ab - 1|
function weighted_psi_symmetric(profile::Vector{<:Dict{Symbol,Int}})
    # Collect all candidates from the first ranking (assuming all share the same keys)
    candidates = collect(keys(profile[1]))
    # All unordered pairs of distinct candidates
    candidate_pairs = collect(combinations(candidates, 2))

    measure = 0.0
    for (c1, c2) in candidate_pairs
        p_ab = ptilde(c1, c2, profile)          # Weighted proportion c1 ≻ c2
        measure += (1.0 - abs(2.0 * p_ab - 1.0)) # 1 - |2 p_ab - 1|
    end

    return measure / length(candidate_pairs)     # Normalize by number of pairs
end



function get_paired_rankings_and_proportions(profile)
    
    tupled = map(x->collect(x) .|> Tuple , profile)
    sorted_tupled = map(x->sort(x; by = x -> x[2]), tupled)
    rankings = [Tuple([string(x[1]) for x in ranking_origin]) for ranking_origin in sorted_tupled]
    
    unique_rankings = unique(rankings)
    
    test_paired, _ = find_reversal_pairs(unique_rankings)
    proportion_test_rankings = proportionmap(rankings)

    return(test_paired, proportion_test_rankings)
end


    
function calc_total_reversal_component(profile)
         paired, proportion_rankings = get_paired_rankings_and_proportions(profile)
          return isempty(paired) ? 0.0 : calc_total_reversal_component(paired, proportion_rankings)
end
    

function calc_reversal_HHI(profile)
    paired, proportion_rankings =  get_paired_rankings_and_proportions(profile)
    return isempty(paired) ? 0.0 : calc_reversal_HHI(paired, proportion_rankings)
end

function  fast_reversal_geometric(profile)
    paired, proportion_rankings =  get_paired_rankings_and_proportions(profile)
    geometric_reversal = fast_reversal_geometric(paired, proportion_rankings)
    return(geometric_reversal)
end



function consensus_for_group(subdf)
    # subdf is a subset of df for a particular (religion, race) group
    # Extract the profiles from this sub-dataframe
    group_profiles = collect(subdf.profile)  # Vector{Dict{Symbol, Int}}
    
    # Call your existing function
    consensus_back_from_permallows, consensus_dict = get_consensus_ranking(group_profiles)
    
    # Return as a NamedTuple so DataFrames can handle it easily
    return (consensus_ranking = consensus_dict)
end



"""
    kendall_tau_dict(r1::Dict{T,Int}, r2::Dict{T,Int}) where T

Given two ranking dictionaries (mapping candidate → rank, with lower numbers = better),
returns the Kendall tau distance (number of discordant candidate pairs).

Two candidates (a, b) are discordant if
    (r1[a] < r1[b]) != (r2[a] < r2[b]).
"""
function kendall_tau_dict(r1::Dict{T,Int}, r2::Dict{T,Int}) where T
    keys_ = collect(keys(r1))
    n = length(keys_)
    d = 0
    for i in 1:(n-1)
        for j in i+1:n
            a = keys_[i]
            b = keys_[j]
            if ( (r1[a] < r1[b]) != (r2[a] < r2[b]) )
                d += 1
            end
        end
    end
    return d
end

"""
    kendall_tau_perm(p1::AbstractVector{T}, p2::AbstractVector{T}) where T

Given two permutations of the same candidate set (each an ordered collection of candidates),
computes the Kendall tau distance (number of discordant pairs).

This function constructs an index mapping for `p2`, then maps `p1` accordingly,
and counts the number of inversions in the resulting array.
"""
function kendall_tau_perm(p1, p2::AbstractVector{T}) where T
    # Build a mapping: candidate -> position in p2.
    pos = Dict{T,Int}()
    for (i, cand) in enumerate(p2)
         pos[cand] = i
    end
    # Map each candidate in p1 to its position in p2.
    mapped = [ pos[cand] for cand in p1 ]
    
    # Count inversions in the mapped array.
    d = 0
    n = length(mapped)
    for i in 1:(n-1)
        for j in i+1:n
            if mapped[i] > mapped[j]
                d += 1
            end
        end
    end
    return d
end


"""
    average_normalized_distance(profile, consensus)

Given a profile (a Vector of Dict{Symbol,Int}) and a consensus ranking `consensus` (also a Dict),
computes the average normalized Kendall tau distance between each ranking in the profile and `consensus`.

The normalization factor is binomial(m,2) where m is the number of alternatives.
"""
function average_normalized_distance(profile, consensus)
    n = length(profile)
    m = length(consensus)
    norm_factor = binomial(m, 2)
    
    dtau = kendall_tau_dict

    total_distance = sum(dtau(ranking, consensus) for ranking in profile)
    return total_distance / (n * norm_factor)
end

# TODO: scrape this subdf bullshit 
# it should just be a damn dictionary 


function group_avg_distance(subdf)
    # Convert subdf.profile (a SubArray) to a plain Vector of dictionaries.
    group_profiles = collect(subdf.profile)
    # Compute the consensus ranking for this group.
    _, consensus_dict = get_consensus_ranking(group_profiles)
    # Compute the average normalized distance.
    avg_dist = average_normalized_distance(group_profiles, consensus_dict)
    # Compute group coherence.
    group_coherence = 1.0 - avg_dist
    return (avg_distance = avg_dist, group_coherence = group_coherence)
end


function weighted_coherence(results_distance::DataFrame, proportion_map::Dict, key)
    total = sum(row.group_coherence * proportion_map[row[key]] for row in eachrow(results_distance))

    # Cstar = (2*total - 1)
    # 3) floor at zero (any tiny numerical dips below 0 become exactly zero)
    #return max(Cstar, 0.0)
    return total
end



# TODO: make this work 
# TODO: be certain this is actually implementing what I thought of 

function pairwise_group_divergence(profile_i, consensus_j, m::Int)
    n_i = length(profile_i)
    norm_factor = binomial(m, 2)
    dtau = kendall_tau_dict
    # Compute the average normalized distance:
    avg_dist = sum(dtau(r, consensus_j) for r in profile_i) / (n_i * norm_factor)
    return avg_dist
end



# TODO: this seems to damn close to 0. Why? 
"""
    overall_divergence(group_profiles, consensus_map)

Given:
  - group_profiles: Dict{T,Vector{Dict{Symbol,Int}}} mapping group id to its profile.
  - consensus_map: Dict{T,Dict{Symbol,Int}} mapping group id to its consensus ranking.
  
Assumes each ranking involves m alternatives (inferred from any consensus ranking).
Computes the divergence measure

   D = (1/((k-1))) * sum_{i ≠ j} ( (n_i/n) * AvgDist(G_i, ρ_j) )

where AvgDist(G_i, ρ_j) is computed using pairwise_divergence.
"""
function overall_divergence(group_profiles, 
                            consensus_map)
    groups = keys(group_profiles)
    k = length(groups)
    # total number of rankings across all groups
    n = sum(length(profile) for profile in values(group_profiles))
    # infer m from any consensus ranking (assumes all have same number of alternatives)
    m = length(first(values(consensus_map)))
    
    total = 0.0
    for i in groups
        n_i = length(group_profiles[i])
        for j in groups
            if i != j
                # divergence from group i to consensus ranking of group j
                d_ij = pairwise_group_divergence(group_profiles[i], consensus_map[j], m)
               # println("Divergence from group $i to consensus of group $j: $d_ij")
               
               #println("n_i: $n_i, n: $n")
                total += (n_i / n) * d_ij
            end
        end
    end
    #println(total)
    D = total / ((k - 1))
    return D
end


function overall_divergences(grouped_consensus, whole_df, key)
    k = nrow(grouped_consensus)
    groups_profiles = Dict(grouped_consensus[i,key] =>
    map(x->x.profile, Base.filter(x-> x[key] == grouped_consensus[i,key], eachrow(whole_df)))
    for i in 1:k) 
    consensus_map = Dict(i[key] => i.x1 for i in eachrow(grouped_consensus))

    D = overall_divergence(groups_profiles, consensus_map)

return D
end


function compute_coherence_and_divergence(df::DataFrame, key::Symbol)
    # 1) Group the DataFrame by `key`
    grouped_df = groupby(df, key)

    # 2) Compute average distance for each group
    results_distance = combine(grouped_df) do subdf
        group_avg_distance(subdf)  # returns (avg_distance = ..., group_coherence = ...)
    end

    # 3) Compute group proportions, e.g. proportion of each group in the entire df
    group_proportions = proportionmap(df[!, key])  # your user-defined function

    # 4) Compute weighted coherence (sum of group_coherence * group proportion)
    #    or if your code uses the `avg_distance` column, adapt accordingly:
    coherence = weighted_coherence(results_distance, group_proportions, key)

    # 5) Compute the consensus ranking for each group
    #    (assumes `consensus_for_group(subdf)` returns (consensus_ranking = ...))
    grouped_consensus = combine(grouped_df) do subdf
        consensus_for_group(subdf)
    end

    # 6) Compute the overall divergence measure
    divergence = overall_divergences(grouped_consensus, df, key)

    return coherence, divergence
end


function apply_measure_to_bts(bts, measure)
    Dict(x => map(y->measure(y.profile), bts[x]) for x in keys(bts))
end


#= function apply_measure_to_encoded_bts(bts, measure)
    Dict(x => map(y->measure(decode_profile_vec(y.profile,bts[x][:code_to_r])), bts[x][:dfs]) for x in keys(bts))
end
 =#


function apply_measure_to_encoded_bts(bts, measure)
    out = Dict{Symbol,Vector{Float64}}()

    for (variant, enc) in pairs(bts)
        dfs        = enc[:dfs]                # Vector{DataFrame}
        code_to_r  = enc[:code_to_r]

        vals = Vector{Float64}(undef, length(dfs))

        for (i, df) in enumerate(dfs)
            decoded = decode_profile_vec(df.profile, code_to_r)
            vals[i] = measure(decoded)
            empty!(decoded)        # free Vector{Dict} payload
        end

        out[variant] = vals
    end
    return out
end 

#@inline _decode(code::Int, code_to_r)::Dict{Symbol,Int} = code_to_r[code]

#= measures = [Ψ, weighted_psi_symmetric,psi_we, psi_wb, calc_total_reversal_component,
                     calc_reversal_HHI, fast_reversal_geometric]
 =#
function apply_all_measures_to_bts(bts; measures = [Ψ,  calc_total_reversal_component,
                     calc_reversal_HHI, fast_reversal_geometric])
        Dict(nameof(measure) => apply_measure_to_bts(bts, measure) for measure in measures)
end

function apply_all_measures_to_encoded_bts(bts; measures = [Ψ, calc_total_reversal_component,
                     calc_reversal_HHI, fast_reversal_geometric] )
        Dict(nameof(measure) => apply_measure_to_encoded_bts(bts, measure) for measure in measures)
end


#= 
function dict_to_long(d)
    rows = NamedTuple[]  # growable vector of NamedTuples
    for (measure, imps) in d
        for (variant, vec) in imps
            for v in vec
                push!(rows, (measure = string(measure),
                             variant = string(variant),
                             value   = float(v)))
            end
        end
    end
    return DataFrame(rows)
end


 =#



function compute_group_metrics(df::DataFrame, demo)
    g = groupby(df, demo)
    results_distance = combine(g) do subdf
        group_avg_distance(subdf)
    end
 
    prop = proportionmap(df[!, demo])
    C = weighted_coherence(results_distance, prop, demo)

    consensus = combine(g) do subdf
        consensus_for_group(subdf)
    end

  
    D = overall_divergences(consensus, df, demo)

    return C, D
end

"""
    compute_demographic_metrics(df::DataFrame,
                                demo_map::Dict{Symbol,String})

Given `df` with `:profile` plus any number of demographics, `demo_map`
maps each column Symbol to the output label.  Returns a DataFrame
with rows (Demographic, C, D).
"""
function compute_demographic_metrics(df::DataFrame,
                                     demo_map::Dict{Symbol,String})
    result = DataFrame(Demographic=String[], C=Float64[], D=Float64[])
    for (col_sym, label) in demo_map
        @assert hasproperty(df, col_sym) "no column $col_sym in df"
        C, D = compute_group_metrics(df, col_sym)
        push!(result, (label, C, D))
    end
    return result
end





function group_consensus_and_divergence(df::DataFrame, demo::Symbol)
    # 1) group and collect keys in order
    gdf        = groupby(df, demo)
    group_keys = [subdf[1, demo] for subdf in gdf]

    # 2) proportion of each group
    prop_map = proportionmap(df[!, demo])

    # 3) build profiles and consensus maps
    group_profiles = Dict(
        k => collect(subdf.profile)
        for (k, subdf) in zip(group_keys, gdf)
    )
    consensus_map = Dict(
        k => get_consensus_ranking(group_profiles[k])[2]
        for k in group_keys
    )

    # 4) build the consensus_df with avg_distance and proportion
    consensus_df = DataFrame()
    consensus_df[!, demo]               = group_keys
    consensus_df[!, :consensus_ranking] = [consensus_map[k] for k in group_keys]
    consensus_df[!, :avg_distance]      = [
        group_avg_distance(subdf).avg_distance
        for subdf in gdf
    ]
    consensus_df[!, :proportion]        = [prop_map[k] for k in group_keys]

    # 5) compute pairwise divergences
    m     = length(first(values(consensus_map)))
    klen  = length(group_keys)
    M     = zeros(Float64, klen, klen)
    for i in 1:klen, j in 1:klen
        M[i,j] = i == j ? 0.0 :
                 pairwise_group_divergence(
                   group_profiles[group_keys[i]],
                   consensus_map[group_keys[j]],
                   m
                 )
    end

    # 6) build the divergence_df with proportion
    col_syms     = Symbol.(string.(group_keys))
    columns_dict = Dict(col_syms[j] => M[:,j] for j in 1:klen)
    divergence_df = DataFrame(columns_dict)
    divergence_df[!, demo]      = group_keys
    divergence_df[!, :proportion] = [prop_map[k] for k in group_keys]
    select!(divergence_df, [demo, :proportion, col_syms...])

    return consensus_df, divergence_df
end







function bootstrap_demographic_metrics(bt_profiles::Dict,
    demo_map::Dict{Symbol,String})

out = Dict{Symbol, Dict{Symbol, Vector{DataFrame}}}()

for (variant, reps) in bt_profiles
C_list = Vector{DataFrame}()
D_list = Vector{DataFrame}()

for df in reps
mdf = compute_demographic_metrics(df, demo_map)   # Demographic | C | D

push!(C_list,  mdf[:, [:Demographic, :C]])
push!(D_list,  mdf[:, [:Demographic, :D]])
end

out[variant] = Dict(:C => C_list,
:D => D_list)
end

return out
end


function bootstrap_group_metrics(bt_profiles, demo)
    result = Dict{Symbol, Dict{Symbol, Vector{Float64}}}()

    for (variant, reps) in bt_profiles
        Cvals = Float64[]
        Dvals = Float64[]

        for df in reps
            C, D = compute_group_metrics(df, demo)
            push!(Cvals, C)
            push!(Dvals, D)
        end

        result[variant] = Dict(:C => Cvals,
                               :D => Dvals)
    end

    return result
end
#= 

function encoded_bootstrap_group_metrics(bts::Dict, demo::Symbol)
result = Dict{Symbol, Dict{Symbol, Vector{Float64}}}()
for (variant, enc) in pairs(bts)
    dfs       = enc[:dfs]          # Vector{DataFrame}
    code_to_r = enc[:code_to_r]

    Cvals = Float64[]
    Dvals = Float64[]

    for df in dfs
        codes = df.profile    
        println("codes: ", codes)
        #println("codes: ", codes)
        df.profile = decode_profile_vec(codes, code_to_r)

        C, D = compute_group_metrics(df, demo)

        push!(Cvals, C)
        push!(Dvals, D)

        df.profile = codes                     # restore & free decoded
    end

    result[variant] = Dict(:C => Cvals, :D => Dvals)
end
return result
end =#

function encoded_bootstrap_group_metrics(bts::Dict, demo::Symbol)
    out = Dict{Symbol, Dict{Symbol, Vector{Float64}}}()

    for (variant, enc) in pairs(bts)
        c2r   = enc[:code_to_r]
        Cvals, Dvals = Float64[], Float64[]

        for df in enc[:dfs]
            decoded = decode_profile_vec(df.profile, c2r)

            # build *new* 2-column frame: profile + the grouping column
            tmp = DataFrame(profile = decoded)
            tmp[!, demo] = df[!, demo]        # copy just this one column

            C, D = compute_group_metrics(tmp, demo)
            push!(Cvals, C)
            push!(Dvals, D)
        end

        out[variant] = Dict(:C => Cvals, :D => Dvals)
    end
    return out
end 
# works for any Dict‐like ranking

"""
    add_G!(stats)

Given the nested dictionary returned by `bootstrap_group_metrics`

    Dict(
        :zero   => Dict(:C => Vector, :D => Vector),
        :random => Dict(...),
        :mice   => Dict(...)
    )

compute the element-wise geometric mean  G = √(C · D)  for every variant and
store it in the inner dictionary under the key **`:G`**.  
The function mutates `stats` and also returns it for convenience.
"""
function add_G!(stats::Dict{Symbol, Dict{Symbol, Vector{Float64}}})
    for sub in values(stats)
        C = sub[:C]
        D = sub[:D]
        @assert length(C) == length(D) "C and D vectors must be the same length"
        sub[:G] = sqrt.(C .* D)
    end
    return stats
end









# *signed* margin, normalised by electorate size
margin(c₁, c₂, profile) = (nab(c₁, c₂, profile) - nab(c₂, c₁, profile)) / length(profile)

function canonical_pair(a, b)
    a < b ? (a, b) : (b, a)
end

function pairwise_margins(profile)
    cs = sort(collect(keys(profile[1])))             # deterministic order
    res = Dict{Tuple{Symbol,Symbol},Float64}()

    for v in combinations(cs, 2)                     # v is Vector{Symbol}
        a, b       = v                              # unpack once
        res[(a,b)] = margin(a, b, profile)   # key is Tuple
    end
    return res
end


function margins_for_rep(df, code_to_r)
    prof = decode_profile_vec(df.profile, code_to_r)  # existing helper of yours
    return pairwise_margins(prof)
end


function margins_over_bootstrap(bt_variant::Dict)
    dfs        = bt_variant[:dfs]             # Vector{DataFrame}
    code_to_r  = bt_variant[:code_to_r]
    return map(df -> margins_for_rep(df, code_to_r), dfs)
end


function summarize_margins(margin_dicts)
    #println(margin_dicts)
    pairs = keys(first(margin_dicts))
    stats = Dict{Tuple{Symbol,Symbol},Tuple{Float64,Float64}}()
    for p in pairs
        vals = [d[p] for d in margin_dicts]
        stats[p] = (mean(vals), std(vals))
    end
    return stats
end


function margins_dataframe(stats)
    pairs = collect(keys(stats))
    DataFrame(; cand1 = first.(pairs),
                 cand2 = last.(pairs),
                 mean_margin = first.(values(stats)),
                 sd_margin   = last.(values(stats)))
end


"""
    margin_stats(bt_top) → Dict{Symbol,DataFrame}

For every imputation variant (`:zero`, `:random`, `:mice`) contained in
`bt_top`, returns a DataFrame with columns
`cand1, cand2, mean_margin, sd_margin`.
"""
function margin_stats(bt_top::Dict{Symbol,Dict})
    out = Dict{Symbol,DataFrame}()
    for (var, bt_variant) in pairs(bt_top)
        mdicts = margins_over_bootstrap(bt_variant)
        stats  = summarize_margins(mdicts)
        out[var] = margins_dataframe(stats)
    end
    return out
end



##############################################################################
#  Extreme–weighted Can-Özkes-Storcken index Ψ_we
#  – linear shape  g(k) = 1 − (k−1)/(κ−1)   (β = 1)
#  – product rule  w_i(a,b) = g(e_i(a))·g(e_i(b))
#  – returns a value in [0,1]
##############################################################################

extremeness(rank::Int, m::Int) = min(rank, m + 1 - rank)

"""
    g_linear_no_cut(k, κ)
Linear decay (★): g(k) = 1 - (k-1)/κ, always positive.
"""
@generated function g_linear_no_cut(k::Int, κ::Int)
    :( 1.0 - (k - 1) / κ )
end

function w_extreme(c1::Symbol, c2::Symbol, ranking::Dict{Symbol,Int})
    m = length(ranking)
    κ = fld(m, 2)
    g1 = g_linear_no_cut(extremeness(ranking[c1], m), κ)
    g2 = g_linear_no_cut(extremeness(ranking[c2], m), κ)
    return g1 * g2                     # product rule
end

function ptilde_extreme(c1::Symbol, c2::Symbol,
                        profile::Vector{<:Dict{Symbol,Int}})
    num = 0.0;  den = 0.0
    for ranking in profile
        w = w_extreme(c1, c2, ranking)
        den += w
        ranking[c1] < ranking[c2] && (num += w)
    end
    return num / den                   # den > 0 for every pair now
end

function psi_we(profile::Vector{<:Dict{Symbol,Int}})
    @assert !isempty(profile)
    cand_pairs = combinations(collect(keys(profile[1])), 2)
    s = 0.0
    for (c1, c2) in cand_pairs
        p̃ = ptilde_extreme(c1, c2, profile)
        s += 1.0 - abs(2p̃ - 1)
    end
    return s / length(cand_pairs)
end






##############################################################################
#  Bottom-weighted Can-Özkes-Storcken index  Ψ_wb
#  • linear “bottom intensity”  h(k) = (k − 1)/(m − 1)   (top = 0, bottom = 1)
#  • product rule  w_i(a,b) = h(rank_i(a)) · h(rank_i(b))
#  • pairs with total weight 0 are ignored in the average
##############################################################################

using Combinatorics            # for `combinations`

# ---------- 1. bottom intensity for a single alternative --------------------

"""
    h_bottom(rank, m) -> Float64

Linear weight favouring the bottom of the ballot:

    h(k) = (k - 1)/(m - 1)         for 1 ≤ k ≤ m,  m ≥ 2
    h(1) = 0   (top),   h(m) = 1   (bottom)
"""
h_bottom(rank::Int, m::Int) = (rank - 1) / (m - 1)

# ---------- 2. pair weight (product rule) -----------------------------------

"""
    w_bottom(c1, c2, ranking) -> Float64

Weight that voter `ranking` assigns to the unordered pair {c1,c2}.
`ranking` is a Dict{Symbol,Int} mapping candidate → rank.
"""
function w_bottom(c1::Symbol, c2::Symbol, ranking::Dict{Symbol,Int})
    m  = length(ranking)
    h1 = h_bottom(ranking[c1], m)
    h2 = h_bottom(ranking[c2], m)
    return h1 * h2
end

# ---------- 3. weighted comparison proportion -------------------------------

"""
    ptilde_bottom(c1, c2, profile) -> Union{Nothing, Float64}

Weighted proportion of voters preferring `c1` to `c2`, using the bottom rule.
Returns `nothing` if the total weight for the pair is zero (then the pair
is skipped in the polarisation average).
"""
function ptilde_bottom(c1::Symbol, c2::Symbol,
                       profile::Vector{<:Dict{Symbol,Int}})
    num = 0.0
    den = 0.0
    for ranking in profile
        w = w_bottom(c1, c2, ranking)
        den += w
        ranking[c1] < ranking[c2] && (num += w)
    end
    return den == 0.0 ? nothing : num / den
end

# ---------- 4. Ψ_wb  ---------------------------------------------------------

"""
    psi_wb(profile) -> Float64

Bottom-weighted polarisation index Ψ_wb.
`profile` is a vector of rankings (Dict{Symbol,Int}), one per voter.

* Range: 0 ≤ Ψ_wb ≤ 1
* Ψ_wb = 0 for unanimous ballots
* Ψ_wb = 1 for the equiprobable profile (all m! permutations once)
"""
function psi_wb(profile::Vector{<:Dict{Symbol,Int}})
    @assert !isempty(profile) "profile must contain at least one ranking"

    candidates      = collect(keys(profile[1]))
    candidate_pairs = collect(combinations(candidates, 2))

    sum_terms = 0.0
    counted   = 0
    for (c1, c2) in candidate_pairs
        p̃ = ptilde_bottom(c1, c2, profile)
        p̃ === nothing && continue             # skip weight-zero pair
        sum_terms += 1.0 - abs(2p̃ - 1)
        counted   += 1
    end
    return counted == 0 ? 0.0 : sum_terms / counted
end


