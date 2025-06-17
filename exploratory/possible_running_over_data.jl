# Refactored Modular ESEB Analysis Code
# =====================================

using Revise
using PrefPol

import PrefPol as pp

using RCall
using LaTeXStrings
using PyCall 
using DataFrames, Statistics
using Colors, ColorSchemes




# This works fine


#= elec = pp.load_election_cfg("config/2022.toml")   # full file into memory
test = pp.weighted_bootstrap(elec)



raw  = pp.load_election_data(elec)               # one call, no branches


bc   = pp.make_bootstrap_cfg(elec)

f,fpath = pp.save_bootstrap(elec; overwrite = true)  # save the bootstrap to disk


 =#

f2 = pp.save_all_bootstraps()


f3 = pp.load_all_bootstraps()

f3[2006].data[1].Lula |> unique

# 1) Impute & save everything you already bootstrapped
paths = pp.impute_all_bootstraps()

#= # 2) Only 2010 and 2014, forcing regeneration, and using a
#    custom ‘most_known_candidates’ vector
paths_sel = impute_all_bootstraps(
               years = [2010, 2014],
               overwrite = true,
               most_known_candidates = ["LULA", "BOLSONARO", "CIRO_GOMES"]
           )

# 3) One single bundle you already have in memory
bt2022 = load_all_bootstraps(years = 2022)[2022]
file   = impute_and_save(bt2022) =#


imps = pp.load_all_imputed_bootstraps()




pp.pp_proportions(f3[2018].data[1], [:Fernando_Haddad])

pp.pp_proportions(imps[2018].data[:mice][1], [:Fernando_Haddad])

imps[2022].data[:zero][1]


elec06 = pp.load_election_cfg("config/2006.toml") 

df_06 = pp.load_election_data(elec06)

df_06.eseb8 |> pp.proportionmap


# for every year, scenario, m there should be a single most_known vector 
# I must check that 


most_known2 = pp.compute_candidate_set(f3[2018].data[1],
    candidate_cols  = elec18.candidates,
    m = elec18.max_candidates, force_include = elec18.scenarios[2].candidates)


# for every year, scenario, m there should be a single most_known vector 
# I must check that 

most_knowns = map(d -> pp.compute_candidate_set(d,
    candidate_cols  = elec18.candidates,
    m = elec18.max_candidates, force_include = elec18.scenarios[2].candidates)
,f3[2018].data)


unique(most_knowns)
 
test  = pp.profile_dataframe(imps[2018].data[:zero][1], score_cols = unique(most_knowns)[1],
demo_cols = elec18.demographics)


another_test = pp.generate_profiles_for_year(2006, f3[2006], imps[2006])

another_test["no_forcing"][3][:zero]


profiles = pp.generate_all_profiles(f3, imps)

## imp20 = pp.load_imputed_bootstrap(2022; quiet = true)


#= 
profiles = save_or_load_all_profiles(f3, imps)

# 2. Force regeneration (e.g. after pipeline changes)
profiles = save_or_load_all_profiles(f3, imps; overwrite = true)

# 3. Only build 2018 & 2022
profiles = save_or_load_all_profiles(f3, imps;
                                     years = [2018, 2022])

# 4. Later, in a new session, just load:
profiles = load_all_profiles() =#


# 1) One year, create-or-load:
p2018 = save_or_load_profiles_for_year(2018, f3, imps)

# 2) Entire set (skipping years not bootstrapped):
allp   = pp.save_or_load_all_profiles_per_year(f3, imps)

# 3) Later, in a fresh REPL, just load one:
p2022  = pp.load_profiles_for_year(2022)


p2018 = pp.load_profiles_for_year(2018)

p2022["lula_bolsonaro"]

p2022["lula_bolsonaro"][2]

pp.apply_all_measures_to_bts(p2022["lula_bolsonaro"][2])


pp.apply_all_measures_to_bts(p2022["lula_bolsonaro"][6])

meas2022 = pp.apply_measures_for_year(p2022)


# 2) Batch for all years:
all_meas = pp.save_or_load_all_measures_per_year(allp)

# 3) Later, just load one:
m2018 = pp.load_measures_for_year(2018)


pp.compute_group_metrics(p2022["lula_bolsonaro"][4][:zero][1], :Ideology)


pp.bootstrap_group_metrics(p2022["lula_bolsonaro"][6], :PT)

pp.bootstrap_group_metrics(p2018["main_four"][6], :Religion)



profiles = pp.save_or_load_all_profiles_per_year(f3, imps)

foo = pp.apply_group_metrics_all_years(allp, f3)

all_gm = pp.save_or_load_all_group_metrics_per_year(allp, f3)



# ================ trying to plot now 




f3       = pp.load_all_bootstraps()                             # year ⇒ (data,cfg,path)
imps     = pp.load_all_imputed_bootstraps()                     # year ⇒ (data,cfg,path)
profiles = pp.save_or_load_all_profiles_per_year(f3, imps)      # year ⇒ scenario ⇒ m ⇒ variant ⇒ [DF]
all_meas = pp.save_or_load_all_measures_per_year(profiles)      # year ⇒ scenario ⇒ m ⇒ measure ⇒ variant ⇒ [values]
all_gm   = pp.save_or_load_all_group_metrics_per_year(profiles, f3)  # year ⇒ scenario ⇒ m ⇒ dem ⇒ measure ⇒ variant ⇒ [values]



all_meas[2006]

# 2) Pick the scenario you care about (exact name must match one in each year):

fi4 = pp.plot_scenario_year(2006, "lula_alckmin", f3, all_meas; variant="mice")

fig5 = pp.plot_scenario_year(2006, "lula_alckmin", f3, all_meas; variant="random")

fig6 = pp.plot_scenario_year(2006, "lula_alckmin", f3, all_meas; variant="zero")


fig = pp.plot_scenario_year(2022, "lula_bolsonaro", f3, all_meas; variant="mice")

fig2 = pp.plot_scenario_year(2022, "lula_bolsonaro", f3, all_meas; variant="random")

fig3 = pp.plot_scenario_year(2022, "lula_bolsonaro", f3, all_meas; variant="zero")


fi7 = pp.plot_scenario_year(2018, "main_four", f3, all_meas; variant="mice")

fig8 = pp.plot_scenario_year(2018, "main_four", f3, all_meas; variant="random")

fig9 = pp.plot_scenario_year(2018, "main_four", f3, all_meas; variant="zero")



f3[2006].data[1][1:20, :] 

imps[2006].data[:zero][1][1:20,:]


imps[2006].data[:mice][1][1:20,:] # the imputation doesn't seem incorrect! 


profiles[2006]["lula_alckmin"][6][:zero][1][1:20,:]



profiles[2006]["lula_alckmin"][6][:zero][1][1,:profile] # the linearization seems to be correct too 

#pp.pp_proportions(f3[2006].data[1], f3[2006].cfg.candidates)

# TODO: test with just two! 


all_gm[2022]["lula_bolsonaro"][2][:PT][:mice]


foo = profiles[2022]["lula_bolsonaro"][6][:mice][1]

gdf = pp.groupby(foo, :PT)

group_keys = [subdf[1, :PT] for subdf in gdf]

    
prop_map = pp.proportionmap(foo[!,:PT])




  # 3) build profiles and consensus maps
group_profiles = Dict(
        k => collect(subdf.profile)
        for (k, subdf) in zip(group_keys, gdf)
    )




    
pp.get_consensus_ranking(group_profiles[0.0])

pp.get_consensus_ranking(group_profiles[1.0])

pp.get_consensus_ranking(group_profiles[99.])


pp.proportionmap(group_profiles[0.0])

pp.proportionmap(group_profiles[1.0])

pp.proportionmap(group_profiles[99.])


consensus_map = Dict(
        k => pp.get_consensus_ranking(group_profiles[k])[2]
        for k in group_keys
    )

# until here it works perfectly!
# it seems to vindicate combining 99 with 1 ! (no, the proportions differ a lot!)




consensus_df = pp.DataFrame()

consensus_df[!, :PT]               = group_keys

consensus_df[!, :consensus_ranking] = [consensus_map[k] for k in group_keys]

consensus_df[!, :avg_distance]      = [
        pp.group_avg_distance(subdf).avg_distance
        for subdf in gdf
    ]

consensus_df[!, :proportion]        = [prop_map[k] for k in group_keys]


avd = 1-pp.average_normalized_distance(group_profiles[1.0],consensus_map[1.0])


    # 5) compute pairwise divergences
m     = length(first(values(consensus_map)))

klen  = length(group_keys)
M     = zeros(Float64, klen, klen)
for i in 1:klen, j in 1:klen
        M[i,j] = i == j ? 0.0 :
                 pp.pairwise_group_divergence(
                   group_profiles[group_keys[i]],
                   consensus_map[group_keys[j]],
                   m
                 )
end

# 6) build the divergence_df with proportion
col_syms     = Symbol.(string.(group_keys))
columns_dict = Dict(col_syms[j] => M[:,j] for j in 1:klen)
divergence_df = pp.DataFrame(columns_dict)
divergence_df[!, :PT]      = group_keys
divergence_df[!, :proportion] = [prop_map[k] for k in group_keys]
pp.select!(divergence_df, [:PT, :proportion, col_syms...])

divergence_df






foo2 = profiles[2022]["lula_bolsonaro"][4][:mice][1]


pp.compute_coherence_and_divergence(foo2, :Ideology)



all_gm[2022]["lula_bolsonaro"][4][:Ideology][:mice]