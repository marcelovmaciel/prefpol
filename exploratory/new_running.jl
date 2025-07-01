using Revise

using PrefPol

import PrefPol as pp



# f2 = pp.save_all_bootstraps()




f3       = pp.load_all_bootstraps()                             # year ⇒ (data,cfg,path)



#paths = pp.impute_all_bootstraps()

#imps2018     = pp.load_all_imputed_bootstraps(years=[2018])                     # year ⇒ (data,cfg,path)

imps     = pp.load_all_imputed_bootstraps()                     # year ⇒ (data,cfg,path)

imps


#profiles = pp.save_or_load_all_profiles_per_year(f3, imps)      # year ⇒ scenario ⇒ m ⇒ variant ⇒ [DF]


#profiles2022 = pp.generate_profiles_for_year(2022, f3[2022], imps[2022])  # scenario ⇒ m ⇒ variant ⇒ [DF]

#profiles2022 = pp.save_or_load_profiles_for_year(2022, f3, imps)  # scenario ⇒ m ⇒ variant ⇒ [DF]

#pp.@save "profiles_2022.jld2" profiles2022

#pp.@load "intermediate_data/profiles/profiles_2022.jld2" profiles2022


#profiles2022



profiles2006  = pp.save_or_load_profiles_for_year(2006, f3, imps[2006])  # scenario ⇒ m ⇒ variant ⇒ [DF]


profiles2018  = pp.save_or_load_profiles_for_year(2018, f3, imps)  # scenario ⇒ m ⇒ variant ⇒ [DF]



#all_meas = pp.save_or_load_all_measures_per_year(profiles)      # year ⇒ scenario ⇒ m ⇒ measure ⇒ variant ⇒ [values]


meas2006 = pp.save_or_load_measures_for_year(2006, profiles2006) # scenario ⇒ m ⇒ measure ⇒ variant ⇒ [values]


meas2022 = pp.save_or_load_measures_for_year(2022, profiles2022) # scenario ⇒ m ⇒ measure ⇒ variant ⇒ [values]

#all_gm   = pp.save_or_load_all_group_metrics_per_year(profiles, f3)  # year ⇒ scenario ⇒ m ⇒ dem ⇒ measure ⇒ variant ⇒ [values]


all_gm2006 = pp.save_or_load_group_metrics_for_year(2006, profiles2006, f3[2006])  # scenario ⇒ m ⇒ dem ⇒ measure ⇒ variant ⇒ [values]


all_gm2022 = pp.save_or_load_group_metrics_for_year(2022, profiles2022, f3[2022])  # scenario ⇒ m ⇒ dem ⇒ measure ⇒ variant ⇒ [values]

helper = Dict(2022=>meas2022)

fig = pp.plot_scenario_year(2022, "lula_bolsonaro", f3, helper; variant="mice")



fig2 = pp.plot_scenario_year(2022, "lula_bolsonaro", f3, helper; variant="random")

fig3 = pp.plot_scenario_year(2022, "lula_bolsonaro", f3, helper; variant="zero")


cfg2022 = f3[2022].cfg          # same cfg you used to build `fig`

pp.save_plot(fig,  2022, "lula_bolsonaro", cfg2022; variant = "mice")


pp.save_plot(fig2, 2022, "lula_bolsonaro", cfg2022; variant = "random")

pp.save_plot(fig3, 2022, "lula_bolsonaro", cfg2022; variant = "zero")




 fig4 = pp.plot_group_demographics_lines(Dict(2022=> all_gm2022),
  f3, 2022, "lula_bolsonaro", variants = [:mice], maxcols = 3, clist_size = 120)

pp.save_plot(fig4, 2022, "lula_bolsonaro", cfg2022; variant = "mice")  

#= profiles2022 = nothing 
all_gm2022 = nothing
meas2022 = nothing =#



helper = Dict(2006=>meas2006)

fig = pp.plot_scenario_year(2006, "lula_alckmin", f3, helper; variant="mice")



fig2 = pp.plot_scenario_year(2006, "lula_alckmin", f3, helper; variant="random")

fig3 = pp.plot_scenario_year(2006, "lula_alckmin", f3, helper; variant="zero")



cfg2006 = f3[2006].cfg          # same cfg you used to build `fig`

pp.save_plot(fig,  2006, "lula_alckmin", cfg2006; variant = "mice")


pp.save_plot(fig2, 2006, "lula_alckmin", cfg2006; variant = "random")


pp.save_plot(fig3, 2006, "lula_alckmin", cfg2006; variant = "zero")





fig4 = pp.plot_group_demographics_lines(Dict(2006=> all_gm2006),
  f3, 2006, "lula_alckmin", variants = [:mice], maxcols = 3, clist_size = 120)

pp.save_plot(fig4, 2006, "lula_alckmin", cfg2006; variant = "mice") 



fig = pp.plot_scenario_year(2006, "no_forcing", f3, helper; variant="mice")



fig2 = pp.plot_scenario_year(2006, "no_forcing", f3, helper; variant="random")

fig3 = pp.plot_scenario_year(2006, "no_forcing", f3, helper; variant="zero")



cfg2006 = f3[2006].cfg          # same cfg you used to build `fig`

pp.save_plot(fig,  2006, "no_forcing", cfg2006; variant = "mice")


pp.save_plot(fig2, 2006, "no_forcing", cfg2006; variant = "random")


pp.save_plot(fig3, 2006, "no_forcing", cfg2006; variant = "zero")
