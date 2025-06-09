#  ======== Testing with paradigmatic profiles ============

top_linear_order_profile = bt_top[:mice][:dfs][1].profile

linear_order_profile = bt_all_profiles[:mice][1].profile

# ===== Using eseb2022 data ====== 
Ψ(top_linear_order_profile) # top linear order profile is coming from refactor_eseb2022 

weighted_psi_symmetric(top_linear_order_profile)

calc_total_reversal_component(top_linear_order_profile)

calc_reversal_HHI(top_linear_order_profile)

fast_reversal_geometric(top_linear_order_profile)


# === Using all candidates =====
Ψ(linear_order_profile) # linear order profile is coming from refactor_eseb2022

weighted_psi_symmetric(linear_order_profile)

calc_total_reversal_component(linear_order_profile) # this looks fishy


calc_reversal_HHI(linear_order_profile)

fast_reversal_geometric(linear_order_profile)




# ===== Using equiprobable profile =========


function create_equiprobable_profile()
    candidates = [:MARINA_SILVA, :CIRO_GOMES, :BOLSONARO, :LULA, :GERALDO_ALCKMIN]
    # Collect all permutations of these five candidates
    perms = collect(permutations(candidates))
    test_profile = [Dict(zip(perm, 1:5)) for perm in perms]
    return test_profile
end


function create_equiprobable_profile(n::Int)
    # 1) Generate symbolic candidates like :C1, :C2, ..., :Cn
    candidates = [Symbol("C$(i)") for i in 1:n]

    # 2) Collect all permutations of these candidates
    perms = collect(permutations(candidates))

    # 3) Convert each permutation into a Dict(candidate => rank)
    profile = [Dict(zip(perm, 1:n)) for perm in perms]

    return profile
end

equiprobable_profile = create_equiprobable_profile()

Ψ(equiprobable_profile)

weighted_psi_symmetric(equiprobable_profile) # this shows this doesn't deal with the equiprobable issue 

calc_total_reversal_component(equiprobable_profile)

calc_reversal_HHI(equiprobable_profile)

fast_reversal_geometric(equiprobable_profile)

# === Testing with reversed profile ====


reversed_profile =[ 
    Dict(:MARINA_SILVA => 1, :CIRO_GOMES => 2, :BOLSONARO => 3, :LULA => 4, :GERALDO_ALCKMIN => 5), 
    Dict(:MARINA_SILVA => 5, :CIRO_GOMES => 4, :BOLSONARO => 3, :LULA => 2, :GERALDO_ALCKMIN => 1)
]


Ψ(reversed_profile)

weighted_psi_symmetric(reversed_profile) # this shows this doesn't deal with the equiprobable issue 

calc_total_reversal_component(reversed_profile)

calc_reversal_HHI(reversed_profile)

fast_reversal_geometric(reversed_profile)

# == Other possible shapes == 

top_consensus_bottom_polarization =[ 
        Dict(:MARINA_SILVA => 1, :CIRO_GOMES => 2, :BOLSONARO => 3, :LULA => 4, :GERALDO_ALCKMIN => 5), 
    Dict(:MARINA_SILVA => 1, :CIRO_GOMES => 2, :BOLSONARO => 5, :LULA => 4, :GERALDO_ALCKMIN => 3)
]

bottom_consensus_top_polarization =[ 
    Dict(:MARINA_SILVA => 1, :CIRO_GOMES => 2, :BOLSONARO => 3, :LULA => 4, :GERALDO_ALCKMIN => 5), 
    Dict(:MARINA_SILVA => 3, :CIRO_GOMES => 2, :BOLSONARO => 1, :LULA => 4, :GERALDO_ALCKMIN => 5)
]


Ψ(top_consensus_bottom_polarization)

Ψ(bottom_consensus_top_polarization)


weighted_psi_symmetric(top_consensus_bottom_polarization)

weighted_psi_symmetric(bottom_consensus_top_polarization)


calc_total_reversal_component(bottom_consensus_top_polarization)

calc_reversal_HHI(top_consensus_bottom_polarization)

fast_reversal_geometric(top_consensus_bottom_polarization)

foo = [:A, :B, :C, :D, :E, :F, :G, :H, :I, :J, :K, :L, :M, :N, :O, :P]


# Yet another top bottom examples 
top3_consensus_bottom_polarization = [
    Dict(:A => 1 , :B => 2, :C=>3, :D => 4, :E => 5, :F => 6, :G => 7, :H => 8, :I => 9, :J=>10),
    Dict(:A => 1 , :B => 2, :C=>3, :D => 10, :E => 9, :F => 8, :G => 7, :H => 6, :I => 5, :J=>4)
]

Ψ(top3_consensus_bottom_polarization)

weighted_psi_symmetric(top3_consensus_bottom_polarization)

calc_total_reversal_component(top3_consensus_bottom_polarization)

calc_reversal_HHI(top3_consensus_bottom_polarization)

fast_reversal_geometric(top3_consensus_bottom_polarization)



bottom3_consensus_top_polarization = [
    Dict(:A => 1 , :B => 2, :C=>3, :D => 4, :E => 5, :F => 6, :G => 7, :H => 8, :I => 9, :J=>10),
    Dict(:A => 7 , :B => 6, :C=>5, :D => 4, :E => 3, :F => 2, :G => 1, :H => 8, :I => 9, :J=>10)
]



Ψ(bottom3_consensus_top_polarization)

weighted_psi_symmetric(bottom3_consensus_top_polarization)

fast_reversal_geometric(bottom3_consensus_top_polarization)

# This shows that the weighted measure is sensitive to where the polarization is happening 

# others 

test_again = [Dict(
    :MARINA_SILVA      => 5,
    :CIRO_GOMES         => 1,
    :BOLSONARO          => 2,
    :LULA               => 3,
    :GERALDO_ALCKMIN    => 4
),
Dict(
    :MARINA_SILVA      => 1,
    :CIRO_GOMES         => 2,
    :BOLSONARO          => 3,
    :LULA               => 4,
    :GERALDO_ALCKMIN    => 5
),
Dict(
    :MARINA_SILVA      => 5,
    :CIRO_GOMES         => 1,
    :BOLSONARO          => 2,
    :LULA               => 3,
    :GERALDO_ALCKMIN    => 4
)]


Ψ(test_again)

weighted_psi_symmetric(test_again)

calc_total_reversal_component(test_again)

yap = [
    Dict(
        :MARINA_SILVA      => 1,
        :CIRO_GOMES         => 2,
        :BOLSONARO          => 3,
        :LULA               => 4,
        :GERALDO_ALCKMIN    => 5
    ),
    Dict(
        :MARINA_SILVA      => 5,
        :CIRO_GOMES         => 4,
        :BOLSONARO          => 3,
        :LULA               => 2,
        :GERALDO_ALCKMIN    => 1
    ),
    Dict(
        :MARINA_SILVA      => 5,
        :CIRO_GOMES         => 4,
        :BOLSONARO          => 3,
        :LULA               => 1,
        :GERALDO_ALCKMIN    => 2
    )
]


Ψ(yap)
weighted_psi_symmetric(yap)


calc_total_reversal_component(yap)
calc_reversal_HHI(yap)
fast_reversal_geometric(yap)



# Now testing the group-based measure! 


df_for_grouping = DataFrame(
    profile  = top_linear_order_profile,
    religion = scores_df.D10,
    race     = scores_df.D12a
)


profile_grouped_religion = groupby(df_for_grouping, :religion)
profile_grouped_race = groupby(df_for_grouping,:race)



results_distance_race = combine(profile_grouped_race) do subdf
    group_avg_distance(subdf)
end

results_distance_religion =  combine(profile_grouped_religion) do subdf
    group_avg_distance(subdf)
end


race_proportion = proportionmap(df_for_grouping[!,:race])

race_coherence = weighted_coherence(results_distance_race,race_proportion, :race)

religion_proportion = proportionmap(df_for_grouping[!, :religion])



religion_coherence = weighted_coherence(results_distance_religion,religion_proportion, :religion)



consensus_grouped_race =  combine(profile_grouped_race) do subdf
    consensus_for_group(subdf)
end

consensus_grouped_religion =  combine(profile_grouped_religion) do subdf
    consensus_for_group(subdf)
end


C_religion = weighted_coherence(results_distance_religion,religion_proportion, :religion)

D_religion = overall_divergences(consensus_grouped_religion, df_for_grouping, :religion)




C_race, D_race = compute_coherence_and_divergence(df_for_grouping, :race)

compute_coherence_and_divergence(df_for_grouping, :religion)



# Now test with obviouslyt polarized profiles. To see if 
# the measure is sensitive to the polarization of the profile. After all. 

polarized_groups = Dict(
    :a=> [Dict(:MARINA_SILVA => 1, :CIRO_GOMES => 2, :BOLSONARO => 3, :LULA => 4, :GERALDO_ALCKMIN => 5),
    Dict(:MARINA_SILVA => 1, :CIRO_GOMES => 2, :BOLSONARO => 3, :LULA => 4, :GERALDO_ALCKMIN => 5),
    Dict(:MARINA_SILVA => 1, :CIRO_GOMES => 2, :BOLSONARO => 3, :LULA => 4, :GERALDO_ALCKMIN => 5),
    Dict(:MARINA_SILVA => 1, :CIRO_GOMES => 2, :BOLSONARO => 3, :LULA => 4, :GERALDO_ALCKMIN => 5)], 
    :b => [Dict(:MARINA_SILVA => 5, :CIRO_GOMES => 4, :BOLSONARO => 3, :LULA => 2, :GERALDO_ALCKMIN => 1),
    Dict(:MARINA_SILVA => 5, :CIRO_GOMES => 4, :BOLSONARO => 3, :LULA => 2, :GERALDO_ALCKMIN => 1),
    Dict(:MARINA_SILVA => 5, :CIRO_GOMES => 4, :BOLSONARO => 3, :LULA => 2, :GERALDO_ALCKMIN => 1),
    Dict(:MARINA_SILVA => 5, :CIRO_GOMES => 4, :BOLSONARO => 3, :LULA => 2, :GERALDO_ALCKMIN => 1)]
)

consensus_map = Dict(:a=> Dict(:MARINA_SILVA => 1, :CIRO_GOMES => 2, :BOLSONARO => 3, :LULA => 4, :GERALDO_ALCKMIN => 5),
:b => Dict(:MARINA_SILVA => 5, :CIRO_GOMES => 4, :BOLSONARO => 3, :LULA => 2, :GERALDO_ALCKMIN => 1))



overall_divergence(polarized_groups, consensus_map)


# simplifyign the process 

compute_coherence_and_divergence(df_for_grouping, :race )



# -- TODO: create one with self-identified ideological positioining! 
# -- TODO: run these with bootstrap to get some confidence intervals.




# -- Integer profile from example 

profile = [
    # 7 voters with A>B>C
    fill(Dict(:A=>1, :B=>2, :C=>3), 7)...,
    
    # 5 voters with A>C>B
    fill(Dict(:A=>1, :C=>2, :B=>3), 5)...,
    
    # 1 voter with C>A>B
    fill(Dict(:C=>1, :A=>2, :B=>3), 1)...,
    
    # 7 voters with C>B>A
    fill(Dict(:C=>1, :B=>2, :A=>3), 7)...,
    
    # 3 voters with B>C>A
    fill(Dict(:B=>1, :C=>2, :A=>3), 3)...,
    
    # 6 voters with B>A>C
    fill(Dict(:B=>1, :A=>2, :C=>3), 6)...,
]


Ψ(profile)
weighted_psi_symmetric(profile)


equi3 = create_equiprobable_profile(3)

equi4 = create_equiprobable_profile(4)

equi7 = create_equiprobable_profile(7)


equi4


Ψ(equi4)

weighted_psi_symmetric(equi4)

(calc_total_reversal_component(equi4)*calc_reversal_HHI(equi4))

calc_reversal_HHI(equi4)

fast_reversal_geometric(equi4)



equi7

factorial(7)

Ψ(equi7)

weighted_psi_symmetric(equi7)

calc_total_reversal_component(equi7)

calc_reversal_HHI(equi7)

calc_reversal_HHI(equi3)




# Unanimity profile generator: returns a profile where every voter has the identical ranking.
function c_profile(m::Int)
    # Generate candidate symbols: :A, :B, :C, … for m alternatives.
    candidates = [Symbol(Char(65 + i - 1)) for i in 1:m]
    # Create the unanimity ranking: candidate in position i gets rank i.
    ranking = Dict(candidates[i] => i for i in 1:m)
    n_voters = 10  # arbitrary number of voters
    # Return a profile: a Vector with 100 identical rankings.
    return fill(ranking, n_voters)
end

# Reversal profile generator: returns a profile where half the voters have a ranking rₛ
# and half have its reversal. (If n_voters is odd, one extra copy goes to the reversed ranking.)
function i_profile(m::Int)
    # Generate candidate symbols as before.
    candidates = [Symbol(Char(65 + i - 1)) for i in 1:m]
    # Standard ranking rₛ: the candidates appear in the given order (position i gets rank i)
    ranking = Dict(candidates[i] => i for i in 1:m)
    # Reversed ranking: reverse the candidate order and assign ranks 1,...,m accordingly.
    ranking_rev = Dict(reverse(candidates)[i] => i for i in 1:m)
    n_voters = 10  # arbitrary number of voters
    n1 = div(n_voters, 2)
    n2 = n_voters - n1
    # Concatenate n1 copies of ranking and n2 copies of ranking_rev.
    return vcat(fill(ranking, n1), fill(ranking_rev, n2))
end

# Example usage:

u3 = c_profile(3)   # For m = 3 (i.e. candidates A, B, C)
i3   = i_profile(3)    # For m = 3 (i.e. candidates A, B, C)
e3 = create_equiprobable_profile(3)

# Map profile names to their values:
profiles = Dict("u3" => u3, "i3" => i3, "e3" => e3)

# List the functions in an array
functions = [Ψ, weighted_psi_symmetric, calc_total_reversal_component, calc_reversal_HHI, fast_reversal_geometric]

# Loop over functions and profiles, printing the function name, profile name, and result.
for f in functions
    # Convert the function to its name string:
    fname = string(nameof(f))
    for (pname, profile) in pairs(profiles)
        result = f(profile)
        println("$(fname)($(pname)) = $(result)")
    end
end



function print_profiles_results(m::Int)
    # Generate the profiles for m alternatives.
    u = c_profile(m)                  # Unanimity profile for m alternatives
    i = i_profile(m)                  # Reversal profile for m alternatives
    e = create_equiprobable_profile(m)  # Equiprobable profile for m alternatives

    # Map profile names to their values.
    # Here we include m in the name (e.g. "u3" for m = 3) for clarity.
    profiles = Dict("u$(m)" => u, "i$(m)" => i, "e$(m)" => e)

   #=  println("u= ", proportionmap(u))
    println("i= ", proportionmap(i))
    println("e= ", proportionmap(e))    =# 

    # List the functions we want to test.
    functions = [Ψ, weighted_psi_symmetric, calc_total_reversal_component,
                 calc_reversal_HHI, fast_reversal_geometric]

    # Loop over functions and profiles, printing the function name, profile name, and result.
    for f in functions
        # Retrieve the function name as a string.
        fname = string(nameof(f))
        for (pname, profile) in pairs(profiles)
            result = f(profile)
            println("$(fname)($(pname)) = $(result)")
        end
    end
    return nothing
end



print_profiles_results(3)


print_profiles_results(4)


print_profiles_results(5)

print_profiles_results(6)

print_profiles_results(7)

print_profiles_results(8)

weighted_psi_symmetric(i_profile(500))

# beginning profile

ranking_to_dict(r::String) = Dict(Symbol(r[i]) => i for i in eachindex(r))

profile = vcat(
    fill(ranking_to_dict("ABCD"), 15),
    fill(ranking_to_dict("DCBA"), 15),
    fill(ranking_to_dict("DABC"), 10),
    fill(ranking_to_dict("ADCB"), 10),
    fill(ranking_to_dict("DBCA"), 10),
    fill(ranking_to_dict("ABCD"), 10),
    fill(ranking_to_dict("BCDA"), 10),
    fill(ranking_to_dict("CBAD"), 10)
)

map(f-> f(profile), [Ψ, weighted_psi_symmetric, calc_total_reversal_component,
                 calc_reversal_HHI, fast_reversal_geometric])


using Random 
# more tests for grouped polarization 

# ------------------------------------------------------------
# helper: make a single-swap (Kendall-tau distance = 1) variant
# ------------------------------------------------------------
"""
    single_swap(rng, ranking_dict) -> Dict

Return a copy of `ranking_dict` where **exactly one pair of
alternatives is swapped**.
"""
function single_swap(rng::AbstractRNG, rd::Dict{T,Int}) where T
    items  = collect(keys(rd))
    i, j   = rand(rng, 1:length(items), 2)     # distinct positions
    while j == i
        j = rand(rng, 1:length(items))
    end
    # build new ranking
    new = copy(rd)
    a, b = items[i], items[j]
    new[a], new[b] = new[b], new[a]
    return new
end

# ------------------------------------------------------------
# main generator
# ------------------------------------------------------------
"""
    make_two_group_profile(n_per_group::Int;
                           rng = Random.GLOBAL_RNG,
                           items = ["a","b","c","d"]) -> (profiles, group_ids)

* `n_per_group` – number of voters in **each** group (total = 2 × n_per_group)
* `items`       – ordered list of alternatives (defaults to `"a"`…`"d"`)
* Returns:

  1. `profiles` – `Vector{Dict{String,Int}}` of length `2n`
  2. `groups`   – `Vector{Int}` with `1` or `2` for each profile
"""
function make_two_group_profile(n_per_group::Int;
                                rng   = Random.GLOBAL_RNG,
                                items = ["a","b","c","d"])

    m = length(items)

    # consensus dictionaries
    cons1 = Dict(items[i] => i for i in 1:m)          # a b c d  → 1 2 3 4
    cons2 = Dict(items[i] => m+1-i for i in 1:m)      # d c b a  → 1 2 3 4

    profiles = Vector{Dict{String,Int}}(undef, 2n_per_group)
    groups   = Vector{Int}(undef, 2n_per_group)

    # fill first group
    for k in 1:n_per_group
        profiles[k]      = single_swap(rng, cons1)
        groups[k]        = 1
    end
    # fill second group
    for k in 1:n_per_group
        profiles[n_per_group + k] = single_swap(rng, cons2)
        groups[n_per_group + k]   = 2
    end

    return profiles, groups
end




function make_two_group_profile_noswap(n_per_group::Int;
                                       items = ["a","b","c","d"])

    m      = length(items)
    cons1  = Dict(items[i] => i         for i in 1:m)     # a b c d
    cons2  = Dict(items[i] => m+1-i     for i in 1:m)     # d c b a

    profiles = [
        (copy(cons1) for _ in 1:n_per_group)...,
        (copy(cons2) for _ in 1:n_per_group)...,
    ]
    group_ids = [fill(1, n_per_group); fill(2, n_per_group)]

    return profiles, group_ids
end


# --- one helper: do *two* adjacent swaps -----------------------------------
@inline function two_adjacent_swaps(rng, order::Vector{T}) where T
    m = length(order)

    i = rand(rng, 1:m-1)         # first adjacent pair i,i+1
    j = rand(rng, 1:m-2)         # second pair j,j+1  (ensure j ≠ i)
    j += (j ≥ i)                 # shift if collision

    order[i], order[i+1] = order[i+1], order[i]
    order[j], order[j+1] = order[j+1], order[j]
    return order
end

# ---------------------------------------------------------------------------
#  main generator
# ---------------------------------------------------------------------------
"""
    make_two_group_profile_twoswap(n_per_group::Int;
                                   rng   = Random.GLOBAL_RNG,
                                   items = ["a","b","c","d"])
        → (profiles, group_ids)

* group 1 consensus = `abcd`
* group 2 consensus = `dcba`
* every ranking differs by **two adjacent swaps** from its group’s consensus
"""
function make_two_group_profile_twoswap(n_per_group::Int;
                                        rng   = Random.GLOBAL_RNG,
                                        items = ["a","b","c","d"])

    m      = length(items)
    cons1  = items                    # a b c d
    cons2  = reverse(items)           # d c b a

    profiles = Vector{Dict{String,Int}}(undef, 2n_per_group)
    groups   = Vector{Int}(undef, 2n_per_group)

    # group 1
    for k in 1:n_per_group
        perm = two_adjacent_swaps(rng, copy(cons1))
        profiles[k] = Dict(c => i for (i,c) in enumerate(perm))
        groups[k]   = 1
    end

    # group 2
    for k in 1:n_per_group
        perm = two_adjacent_swaps(rng, copy(cons2))
        idx  = n_per_group + k
        profiles[idx] = Dict(c => i for (i,c) in enumerate(perm))
        groups[idx]   = 2
    end

    return profiles, groups
end

profiles, g = make_two_group_profile_noswap(10)

profiles, g = make_two_group_profile(300; rng = MersenneTwister(42))

profiles, g = make_two_group_profile_twoswap(300)

df_for_grouping = pp.DataFrame(
    profile  = profiles,
    g = g
)




profile_grouped_g = pp.groupby(df_for_grouping, :g)




results_distance_g = pp.combine(profile_grouped_g) do subdf
    println("subdf: ", subdf)
    pp.group_avg_distance(subdf)
end



g_proportion = pp.proportionmap(df_for_grouping.g)

g_coherence = pp.weighted_coherence(results_distance_g,g_proportion, :g)



consensus_grouped_g =  pp.combine(profile_grouped_g) do subdf
    pp.consensus_for_group(subdf)
end



C_g = pp.weighted_coherence(results_distance_g,g_proportion, :g)
D_g = pp.overall_divergences(consensus_grouped_g, df_for_grouping, :g)





# yet another test 

function make_two_group_profile(; n₁::Int = 90, n₂::Int = 90,
                                   p_swap::Float64 = 0.10)

    cons1 = Dict(:a=>1, :b=>2, :c=>3, :d=>4)
    cons2 = Dict(:a=>4, :b=>3, :c=>2, :d=>1)

    profiles = Vector{Dict{Symbol,Int}}(undef, n₁ + n₂)
    groups   = Vector{Int}(undef, n₁ + n₂)

    # group 1
    for i in 1:n₁
        profiles[i] = copy(cons1);  groups[i] = 1
    end

    # helper – one adjacent swap away from cons2
    one_swap_rank() = let r = copy(cons2)
        top2 = sort(collect(keys(r)); by = c -> r[c])[1:2]   # ← fixed
        r[top2[1]], r[top2[2]] = r[top2[2]], r[top2[1]]
        r
    end

    n_swap = round(Int, p_swap * n₂)

    # group 2
    for j in 1:n₂
        idx = n₁ + j
        profiles[idx] = (j ≤ n_swap) ? one_swap_rank() : copy(cons2)
        groups[idx]   = 2
    end

    return profiles, groups
end


function cdg_for_toy_profile(; kwargs...)
profiles, groups = make_two_group_profile(; kwargs...)
df = pp.DataFrame(profile = profiles, g = groups)
C, D = pp.compute_group_metrics(df, :g)   # <— your existing routine
return C, D, sqrt(C*D)
end

C, D, G = cdg_for_toy_profile(n₁ = 5200, n₂ = 1200, p_swap = 0.10)

C, D, G = cdg_for_toy_profile(n₁ = 200, n₂ = 1200, p_swap = 0.30)

C, D, G = cdg_for_toy_profile(n₁ = 200, n₂ = 1200, p_swap = 0.50)

C, D, G = cdg_for_toy_profile(n₁ = 200, n₂ = 1200, p_swap = 0.700)

C, D, G = cdg_for_toy_profile(n₁ = 200, n₂ = 1200, p_swap = 0.90)

C, D, G = cdg_for_toy_profile(n₁ = 200, n₂ = 1200, p_swap = 1.)




"""
    make_two_group_profile_swaps(;
        n₁::Int        = 90,
        n₂::Int        = 90,
        p_swap₁::Float64 = 0.10,   # % of group 1 that is altered
        p_swap₂::Float64 = 0.10,   # % of group 2 that is altered
        swap_dist::Int   = 1,      # number of adjacent swaps (k ≥ 1)
        rng              = Random.GLOBAL_RNG)

Return `(profiles, groups)` where  

* group 1’s consensus is **a b c d**  
* group 2’s consensus is **d c b a**  

and the requested share of each group is exactly `swap_dist`
adjacent-swaps away from its own consensus.
"""
function make_two_group_profile_swaps(; n₁::Int        = 90,
                                         n₂::Int        = 90,
                                         p_swap₁::Float64 = 0.10,
                                         p_swap₂::Float64 = 0.10,
                                         swap_dist::Int   = 1,
                                         rng              = Random.GLOBAL_RNG)

    cons1 = Dict(:a=>1, :b=>2, :c=>3, :d=>4)   # a b c d
    cons2 = Dict(:a=>4, :b=>3, :c=>2, :d=>1)   # d c b a

    # --- helper -----------------------------------------------------------
    function k_adjacent_swaps(r::Dict{Symbol,Int}, k::Int)
        out = copy(r)
        keys_vec = sort(collect(keys(out)); by = c -> out[c])  # [a,b,c,d]
        for _ in 1:k
            i = rand(rng, 1:length(keys_vec)-1)                # pick position i
            c1, c2 = keys_vec[i], keys_vec[i+1]
            out[c1], out[c2] = out[c2], out[c1]                # swap their ranks
        end
        return out
    end
    # ----------------------------------------------------------------------

    profiles = Vector{Dict{Symbol,Int}}(undef, n₁ + n₂)
    groups   = Vector{Int}(undef, n₁ + n₂)

    # group 1
    n_swap₁ = round(Int, p_swap₁ * n₁)
    for i in 1:n₁
        profiles[i] = (i ≤ n_swap₁) ?
                       k_adjacent_swaps(cons1, swap_dist) :
                       copy(cons1)
        groups[i]   = 1
    end

    # group 2
    n_swap₂ = round(Int, p_swap₂ * n₂)
    for j in 1:n₂
        idx = n₁ + j
        profiles[idx] = (j ≤ n_swap₂) ?
                        k_adjacent_swaps(cons2, swap_dist) :
                        copy(cons2)
        groups[idx]   = 2
    end

    return profiles, groups
end


function cdg_for_toy_profile_swaps(; kwargs...)
    profiles, groups = make_two_group_profile_swaps(; kwargs...)
    df = pp.DataFrame(profile = profiles, g = groups)   # pp = your module
    C, D = pp.compute_group_metrics(df, :g)
    return C, D, sqrt(C*D)
end




C, D, G = cdg_for_toy_profile_swaps(n₁ = 1990, n₂ = 1990,
                                    p_swap₁ = 0.001,
                                    p_swap₂ = 0.99,
                                    swap_dist = 2)


# 200 vs. 300 voters, 20 % and 15 % deviating by *two* adjacent swaps
C2, D2, G2 = cdg_for_toy_profile_swaps(n₁ = 1200, n₂ = 1200,
                                       p_swap₁ = 0.70,
                                       p_swap₂ = 0.7,
                                       swap_dist = 4)




function avg_distance_to(target_consensus, profiles)
    m = length(target_consensus)
    norm = binomial(m,2)
    pp.mean(pp.kendall_tau_dict(r, target_consensus) / norm for r in profiles)
end

# build a profile as before
prof, grp = make_two_group_profile_swaps(n₁=1990, n₂=1990,
                                         p_swap₁=0.001,
                                         p_swap₂=0.99,
                                         swap_dist=4)
cons1 = Dict(:a=>1,:b=>2,:c=>3,:d=>4)
cons2 = Dict(:a=>4,:b=>3,:c=>2,:d=>1)

println("avg distance G1 → ρ2 = ",
        avg_distance_to(cons2, prof[grp .== 1]))  # ≈ 1.0 (max)

println("avg distance G2 → ρ1 = ",
        avg_distance_to(cons1, prof[grp .== 2]))  # ≈ 0.67                                       






"""
    toy_groups(m, n1, n2; dist_consensus, p1, p2, swap_dist=1)

Return `profiles, groups`
  • consensuses are `abcd…` and another permutation at Kendall distance
    `dist_consensus`
  • fraction `p1` of group 1, `p2` of group 2` are `swap_dist`-adjacent-swaps
    away from their own consensus
"""
function toy_groups(m, n1, n2;
                    dist_consensus = 1.0,  # choose one of the distances above
                    p1 = 0.0,
                    p2 = 0.0,
                    swap_dist = 1)

    # 1  consensuses ---------------------------------------------------------
    base = Symbol.('a':'a'+m-1)
    ρ1   = Dict(base[i] => i for i in 1:m)

    # pick a second permutation with the desired normalised distance
    perms = collect(pp.permutations(base))
    norm  = pp.binomial(m,2)
    ρ2vec = first(filter(p -> pp.kendall_tau_perm(p, base)/norm ≈ dist_consensus,
                         perms))
    ρ2    = Dict(ρ2vec[i] => i for i in 1:m)

    consensus = (ρ1, ρ2)

    # helper to make an adjacent-swap copy ------------------------------
    function noisy(r, k)
        cands = sort(collect(keys(r)); by = c -> r[c])
        for _ in 1:k
            i = rand(1:m-1)
            c1, c2 = cands[i], cands[i+1]
            r[c1], r[c2] = r[c2], r[c1]
            cands[i], cands[i+1] = cands[i+1], cands[i]    # keep order array in sync
        end
        r
    end

    profiles = Vector{Dict{Symbol,Int}}()
    groups   = Int[]

    # group 1
    push!(profiles, [copy(ρ1) for _ in 1:round(Int,(1-p1)*n1)]...)
    push!(profiles, [noisy(copy(ρ1), swap_dist) for _ in 1:round(Int,p1*n1)]...)
    append!(groups, fill(1, n1))

    # group 2
    push!(profiles, [copy(ρ2) for _ in 1:round(Int,(1-p2)*n2)]...)
    push!(profiles, [noisy(copy(ρ2), swap_dist) for _ in 1:round(Int,p2*n2)]...)
    append!(groups, fill(2, n2))

    shuffle!(profiles); shuffle!(groups)
    return profiles, groups
end        



###############################################################################
# 1. build a synthetic profile
###############################################################################
# • m          = number of alternatives
# • n1, n2     = group sizes
# • p1, p2     = share of each group that is *swap-dist* swaps away
# • dist_consensus = normalised Kendall-τ distance between the two consensuses
#                    pick one of 0, 1/6, 2/3, 5/6, 1 for m = 4 (see table)
# • swap_dist  = 1 (one adjacent swap), 2, 3 … for more disorder
profiles, groups = toy_groups(4,   # m
                              2000, 2000;          # n1, n2
                              dist_consensus = 1,  # dcba vs abcd
                              p1 = 0.05,           # 5 % noisy in group 1
                              p2 = 0.90,           # 90 % noisy in group 2
                              swap_dist = 1)       # one-swap noise

###############################################################################
# 2. wrap into a DataFrame that your existing pipeline understands
###############################################################################
using DataFrames      # if not already in scope
df = DataFrame(profile = profiles, g = groups)

###############################################################################
# 3. get C, D and G with your existing routine
###############################################################################
C, D = pp.compute_group_metrics(df, :g)   # ← your function
G     = sqrt(C * D)

println("C = ", C)
println("D = ", D)
println("G = ", G)              



# 40 % of both groups are two swaps away
profiles, groups = toy_groups(4, 2000, 2000;
                              dist_consensus = 1,
                              p1 = 0.40, p2 = 0.40,
                              swap_dist = 2)

C, D = pp.compute_group_metrics(DataFrame(profile=profiles, g=groups), :g)
println("C = $C   D = $D   G = $(sqrt(C*D))")



profiles, groups = toy_groups(4, 2000, 2000;
                              dist_consensus = 1/6,
                              p1 = 0.90, p2 = 0.90,
                              swap_dist = 3)

C, D = pp.compute_group_metrics(DataFrame(profile=profiles, g=groups), :g)
println("C = $C   D = $D   G = $(sqrt(C*D))")



function report_toy(m; n1, n2, dist_consensus, p1, p2, swap_dist)
    prof, grp = toy_groups(m, n1, n2;
                           dist_consensus = dist_consensus,
                           p1 = p1, p2 = p2,
                           swap_dist = swap_dist)
    C, D = pp.compute_group_metrics(DataFrame(profile = prof, g = grp), :g)
    G    = sqrt(C*D)
    println(rpad("C = $(round(C, digits=3))",18),
            rpad("D = $(round(D, digits=3))",18),
            "G = $(round(G, digits=3))")
end



m  = 4          # number of alternatives
n1 = n2 = 2000  # size of each group

println("low-C  /  low-D")
report_toy(m; n1, n2, dist_consensus = 1/6,
                p1 = 0.90, p2 = 0.90, swap_dist = 1)

println("medium-C  /  high-D")






report_toy(m; n1, n2, dist_consensus = 1,
                p1 = 0.40, p2 = 0.40, swap_dist = 1)

println("high-C  /  medium-D")
report_toy(m; n1, n2, dist_consensus = 1/3,
                p1 = 0.00, p2 = 0.00, swap_dist = 1)


report_toy(m; n1, n2,
           dist_consensus = 1/6,
           p1 = 0.95, p2 = 0.95,
           swap_dist = 1)
# → C ≈ 0.15   D ≈ 0.15                




report_toy(m; n1, n2,
           dist_consensus = 1/3,   # half-way apart
           p1 = 0.90, p2 = 0.90,
           swap_dist = 1)
# → C ≈ 0.20   D ≈ 0.35



report_toy(m; n1, n2,
           dist_consensus = 1,     # abcd vs dcba
           p1 = 0.90, p2 = 0.90,
           swap_dist = 1)
# → C ≈ 0.18   D ≈ 0.82


report_toy(m; n1, n2,
           dist_consensus = 1/6,
           p1 = 0.45, p2 = 0.45,
           swap_dist = 1)
# → C ≈ 0.52   D ≈ 0.18



report_toy(m; n1, n2,
           dist_consensus = 1/3,
           p1 = 0.40, p2 = 0.40,
           swap_dist = 1)
# → C ≈ 0.50   D ≈ 0.50



report_toy(m; n1, n2,
           dist_consensus = 1,
           p1 = 0.40, p2 = 0.40,
           swap_dist = 1)
# → C ≈ 0.51   D ≈ 0.79


report_toy(m; n1, n2,
           dist_consensus = 1/6,
           p1 = 0.00, p2 = 0.00,   # both groups perfectly follow consensus
           swap_dist = 1)
# → C ≈ 0.83   D ≈ 0.17

report_toy(m; n1, n2,
           dist_consensus = 1/3,
           p1 = 0.00, p2 = 0.00,
           swap_dist = 1)
# → C ≈ 0.83   D ≈ 0.33



report_toy(m; n1, n2,
           dist_consensus = 1,
           p1 = 0.00, p2 = 0.00,
           swap_dist = 1)
# → C ≈ 0.83   D ≈ 0.83



prof, grp = toy_groups(4, 2000, 2000;
                       dist_consensus = 1,     # abcd vs dcba
                       p1 = 0.00,              # group-1 perfectly coherent
                       p2 = 0.80,              # 80 % reverse in group-2
                       swap_dist = 0)          # 0 = full reverse

# -----------------------------------------------------------
# 2. compute C, D, G with your own pipeline
# -----------------------------------------------------------
C, D = pp.compute_group_metrics(DataFrame(profile=prof, g=grp), :g)


G     = sqrt(C*D)



profiles, groups = toy_groups(4, 2000, 2000;
                              dist_consensus = 1,    # abcd ↔ dcba
                              p1 = 0.00,             # G₁ = consensus
                              p2 = 1.00,             # all of G₂ = their consensus
                              swap_dist = 0)         # no extra noise


C, D = pp.compute_group_metrics(DataFrame(profile=profiles, g=groups), :g)


G    = sqrt(C*D)




function toy_groups(m::Int, n1::Int, n2::Int;
                    p1::Float64=0.0,
                    p2::Float64=0.0,
                    swap_dist::Int=1)

    # 1) build the two *fixed* consensuses:
    cons1 = Dict(i => i for i in 1:m)
    cons2 = Dict(i => m+1 - i for i in 1:m)

    # helper: apply k random adjacent‐swaps to ranking r
    function noisy(r::Dict{Int,Int}, k)
        s = sort(collect(r); by = x->x[2])   # (candidate,rank) sorted by rank
        perm = [c[1] for c in s]             # a vector of candidate names
        for _ in 1:k
            # pick a random adjacent position and swap it
            i = rand(1:m-1)
            perm[i], perm[i+1] = perm[i+1], perm[i]
        end
        # re‐encode back into a Dict
        return Dict(c => findfirst(==(c), perm) for c in perm)
    end

    total = n1 + n2
    profiles = Vector{Dict{Int,Int}}(undef, total)
    groups   = Vector{Int}(undef, total)

    # 2) fill group 1
    n1swap = round(Int, p1 * n1)
    for i in 1:n1
        if i ≤ n1swap
            profiles[i] = noisy(cons1, swap_dist)
        else
            profiles[i] = copy(cons1)
        end
        groups[i] = 1
    end

    # 3) fill group 2
    n2swap = round(Int, p2 * n2)
    for j in 1:n2
        idx = n1 + j
        if j ≤ n2swap
            profiles[idx] = noisy(cons2, swap_dist)
        else
            profiles[idx] = copy(cons2)
        end
        groups[idx] = 2
    end

    return profiles, groups
end



profiles, groups = toy_groups(4, 2000, 2000;
                              p1 = 0.0,   # nobody in G1 noisy
                              p2 = 0.0,   # nobody in G2 noisy
                              swap_dist = 0)

using DataFrames
df = DataFrame(profile=profiles, g=groups)

C, D = pp.compute_group_metrics(df, :g)
G     = sqrt(C*D)

@show C, D, G



# Make *both* G₁ and G₂ totally scrambled
profiles, groups = toy_groups(4, 2000, 2000;
                             p1 = 1.0,    # all of G₁ noisy
                             p2 = 1.0,    # all of G₂ noisy
                             swap_dist = 10)

df = DataFrame(profile = profiles, g = groups)
C, D = pp.compute_group_metrics(df, :g)
G     = sqrt(C * D)

println("C = $(round(C,digits=3)), D = $(round(D,digits=3)), G = $(round(G,digits=3))")




function uniform_random_groups(m::Int, n1::Int, n2::Int)
    # 1) label the alternatives however you like; here :a, :b, :c, … 
    alts = Symbol.(string.('a':'a'+(m-1)))
    all_perms = collect(pp.permutations(alts))
    P = length(all_perms)

    # 2) sample n1 + n2 times **with replacement** from the m! permutations
    profiles = Vector{Dict{Symbol,Int}}(undef, n1+n2)
    groups   = Vector{Int}(undef, n1+n2)

    for i in 1:(n1+n2)
        π = all_perms[rand(1:P)]
        # turn π = (a,b,c,...) into Dict(a=>1,b=>2,c=>3,…)
        profiles[i] = Dict(π[j] => j for j in 1:m)
        groups[i] = i ≤ n1 ? 1 : 2
    end

    return DataFrame(profile = profiles, g = groups)
end

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# example parameters
m, n1, n2 = 4, 2000, 2000

df = uniform_random_groups(m, n1, n2)

C, D = pp.compute_group_metrics(df, :g)
G = sqrt(C * D)

println("Uniform–random within each group:")
println("  C = ", round(C, digits=3))
println("  D = ", round(D, digits=3))
println("  G = ", round(G, digits=3))



# reverse a consensus dict
function _reverse_consensus(cons::Dict{Symbol,Int})
    ks = sort(collect(keys(cons)); by = k->cons[k])
    m  = length(ks)
    return Dict(k => (m - i + 1) for (i,k) in enumerate(ks))
end

"""
    perfect_internal_polarized_groups(
        cons::Dict{Symbol,Int},
        n1::Int, n2::Int
    ) -> DataFrame

Builds a DataFrame with two groups of size `n1` and `n2`.  
*Within each group*, half the members are `cons` (e.g. `:a=>1,:b=>2,:c=>3,:d=>4`)  
and half are its exact reversal (e.g. `:a=>4,:b=>3,:c=>2,:d=>1`).
"""
function perfect_internal_polarized_groups(
        cons::Dict{Symbol,Int},
        n1::Int, n2::Int
    )
    rev = _reverse_consensus(cons)
    N   = n1 + n2

    profiles = Vector{Dict{Symbol,Int}}(undef, N)
    groups   = Vector{Int}(undef, N)

    # group 1: first n1 entries
    half1 = fld(n1,2)
    for i in 1:n1
        profiles[i] = (i ≤ half1 ? cons : rev)
        groups[i]   = 1
    end

    # group 2: next n2 entries
    half2 = fld(n2,2)
    for j in 1:n2
        idx = n1 + j
        profiles[idx] = (j ≤ half2 ? cons : rev)
        groups[idx]   = 2
    end

    return DataFrame(profile = profiles, g = groups)
end

# ── usage ─────────────────────────────────────────────────────────────

cons_abcd = Dict(:a=>1, :b=>2, :c=>3, :d=>4)

df = perfect_internal_polarized_groups(cons_abcd, 2000, 2000)

C, D = pp.compute_group_metrics(df, :g)
G     = sqrt(C*D)

@show C, D, G

filter(x->x.g ===1, df )

