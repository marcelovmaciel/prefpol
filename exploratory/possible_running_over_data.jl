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




# 2) Pick the scenario you care about (exact name must match one in each year):

fig = pp.plot_scenario_year(2022, "lula_bolsonaro", f3, all_meas; variant="mice")

fig2 = pp.plot_scenario_year(2022, "lula_bolsonaro", f3, all_meas; variant="random")

fig3 = pp.plot_scenario_year(2022, "lula_bolsonaro", f3, all_meas; variant="zero")