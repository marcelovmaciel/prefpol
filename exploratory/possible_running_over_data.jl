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
p2022  = load_profiles_for_year(2022)



# ===============================================
# from here I need to change 
bc   = pp.make_bootstrap_cfg(elec)
raw  = pp.load_election_data(elec)               # one call, no branches
bt   = pp.build_bootstrapped_encoded_variants(bc)  
# ===============================================




sel_bt = pp.build_boootstraped_scenarios(
            years = [2018, 2022],
            which = Dict(2018 => ["lula_bolsonaro", "main_four"], 2022 => ["lula_bolsonaro"])
         )


all_bt = pp.build_boootstraped_scenarios()         



pp.apply_all_measures_to_encoded_bts(all_bt[2006]["no_forcing"].variants)



all_bt[2006]

# 3. compute evesry global measure for every m
globals = pp.globals_over_m(all_bt[2006]["no_forcing"])

all_bt[2006]["no_forcing"].variants