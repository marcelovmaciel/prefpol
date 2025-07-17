using Revise
using PrefPol

import PrefPol as pp



# f2 = pp.save_all_bootstraps()




f3       = pp.load_all_bootstraps()                             # year ⇒ (data,cfg,path)


#index_paths = pp.impute_from_f3(f3; overwrite=false)


#paths = pp.impute_all_bootstraps()

#imps2018     = pp.load_all_imputed_bootstraps(years=[2018])                     # year ⇒ (data,cfg,path)

iy2006 = pp.load_imputed_year(2006)



iy2018 = pp.load_imputed_year(2018)

iy2022 = pp.load_imputed_year(2022)

profiles2006 = pp.generate_profiles_for_year_streamed_from_index(
                   2006, f3[2006], iy2006; overwrite = false)


profiles2018 = pp.generate_profiles_for_year_streamed_from_index(
                   2018, f3[2018], iy2018; overwrite = false)

profiles2022 = pp.generate_profiles_for_year_streamed_from_index(
                   2022, f3[2022], iy2022; overwrite = false)                   



measures2006 = pp.save_or_load_measures_for_year(2006, profiles2006;
                                              overwrite = false,   # set true to rebuild
                                              verbose   = true)    # progress / info logs


measures2018 = pp.save_or_load_measures_for_year(2018, profiles2018;
                                              overwrite = false,   # set true to rebuild
                                              verbose   = true)    # progress / info logs


measures2022 = pp.save_or_load_measures_for_year(2022, profiles2022;
                                              overwrite = false,   # set true to rebuild
                                              verbose   = true)    # progress / info logs                                              


#= test1 = profiles2022["lula_bolsonaro"][6][:mice, 1]

cs = pp.metadata(test1, "candidates")            # Vector{Symbol}

perm = test1.profile[1]

pp.perm2dict(perm, cs)                 # original Dict ranking


map(x->pp.perm2dict(x, cs), test1.profile) =#


#imps     = pp.load_all_imputed_bootstraps()                     # year ⇒ (data,cfg,path)

imps


#profiles = pp.save_or_load_all_profiles_per_year(f3, imps)      # year ⇒ scenario ⇒ m ⇒ variant ⇒ [DF]


#profiles2022 = pp.generate_profiles_for_year(2022, f3[2022], imps[2022])  # scenario ⇒ m ⇒ variant ⇒ [DF]

#profiles2022 = pp.save_or_load_profiles_for_year(2022, f3, imps)  # scenario ⇒ m ⇒ variant ⇒ [DF]

#pp.@save "profiles_2022.jld2" profiles2022

#pp.@load "intermediate_data/profiles/profiles_2022.jld2" profiles2022


#profiles2022







profiles2006  = pp.save_or_load_profiles_for_year(2006, f3, imps)  # scenario ⇒ m ⇒ variant ⇒ [DF]

imps[2006].data[:mice][1] |> Base.summarysize


profiles2006["lula_alckmin"][6][:mice][1] |> Base.summarysize

profiles2006 |> Base.summarysize



df        = profiles2006["lula_alckmin"][6][:mice][1]


perm      = df.profile[1]                         # SVector{m,UInt8}

cand_syms = pp.metadata(df, "candidates")            # Vector{Symbol}

dict = pp.perm2dict(perm, cand_syms)                 # original Dict ranking

#= 

cands = profiles2006["lula_alckmin"][6][:mice][1][1,:profile]  |> keys |> collect;



df = profiles2006["lula_alckmin"][6][:mice][1]

profiles2006["lula_alckmin"][6][:mice][1][1,:profile] 




pool = pp.compress_rank_column!(df,cands; col = :profile);

profiles2006["lula_alckmin"][6][:mice][1] |> Base.summarysize


sv = df.profile[1];                      # already decoded




pp.perm2dict(sv, cands)   =#



profiles2006




profiles2006["lula_alckmin"][6][:mice]


profiles2006["lula_alckmin"][6][:mice][1]


profiles2006["lula_alckmin"][6][:mice][1][1,:profile] 




profiles2018  = pp.save_or_load_profiles_for_year(2018, f3, imps)  # scenario ⇒ m ⇒ variant ⇒ [DF]




# profiles2022 = pp.save_or_load_profiles_for_year(2022, f3, imps[2022])  # scenario ⇒ m ⇒ variant ⇒ [DF]

pp.@load "intermediate_data/profiles/profiles_2022.jld2" profiles2022  # scenario ⇒ m ⇒ variant ⇒ [DF]


profiles2022

#all_meas = pp.save_or_load_all_measures_per_year(profiles)      # year ⇒ scenario ⇒ m ⇒ measure ⇒ variant ⇒ [values]


meas2006 = pp.save_or_load_measures_for_year(2006, profiles2006) # scenario ⇒ m ⇒ measure ⇒ variant ⇒ [values]


meas2022 = pp.save_or_load_measures_for_year(2022, profiles2022) # scenario ⇒ m ⇒ measure ⇒ variant ⇒ [values]


meas2018 = pp.save_or_load_measures_for_year(2018, profiles2018) # scenario ⇒ m ⇒ measure ⇒ variant ⇒ [values]

#all_gm   = pp.save_or_load_all_group_metrics_per_year(profiles, f3)  # year ⇒ scenario ⇒ m ⇒ dem ⇒ measure ⇒ variant ⇒ [values]


all_gm2006 = pp.save_or_load_group_metrics_for_year(2006, profiles2006, f3[2006])  # scenario ⇒ m ⇒ dem ⇒ measure ⇒ variant ⇒ [values]


all_gm2022 = pp.save_or_load_group_metrics_for_year(2022, profiles2022, f3[2022])  # scenario ⇒ m ⇒ dem ⇒ measure ⇒ variant ⇒ [values]


all_gm2018 = pp.save_or_load_group_metrics_for_year(2018, profiles2018, f3[2018])  # scenario ⇒ m ⇒ dem ⇒ measure ⇒ variant ⇒ [values]

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






helper = Dict(2018=>meas2018)

fig = pp.plot_scenario_year(2018, "main_four", f3, helper; variant="mice")



fig2 = pp.plot_scenario_year(2018, "no_forcing", f3, helper; variant="mice")

fig3 = pp.plot_scenario_year(2018, "lula_bolsonaro", f3, helper; variant="mice")




cfg2018 = f3[2018].cfg          # same cfg you used to build `fig`

pp.save_plot(fig,  2018, "main_four", cfg2018; variant = "mice")


pp.save_plot(fig2, 2018, "no_forcing", cfg2018; variant = "mice")

pp.save_plot(fig3, 2018, "lula_bolsonaro", cfg2018; variant = "mice")





fig4 = pp.plot_group_demographics_lines(Dict(2018=> all_gm2018),
  f3, 2018, "lula_bolsonaro", variants = [:mice], maxcols = 2, clist_size = 60)


pp.save_plot(fig4, 2018, "lula_bolsonaro", cfg2018; variant = "mice") 



fig4 = pp.plot_group_demographics_lines(Dict(2018=> all_gm2018),
  f3, 2018, "lula_bolsonaro", variants = [:zero], maxcols = 2, clist_size = 60)

pp.save_plot(fig4, 2018, "lula_bolsonaro", cfg2018; variant = "zero") 



fig4 = pp.plot_group_demographics_lines(Dict(2018=> all_gm2018),
  f3, 2018, "main_four", variants = [:mice], maxcols = 2, clist_size = 60)


pp.save_plot(fig4, 2018, "main_four", cfg2018; variant = "mice") 



fig = pp.plot_scenario_year(2018, "main_four", f3, helper; variant="zero")



fig2 = pp.plot_scenario_year(2018, "no_forcing", f3, helper; variant="zero")

fig3 = pp.plot_scenario_year(2018, "lula_bolsonaro", f3, helper; variant="zero")




cfg2018 = f3[2018].cfg          # same cfg you used to build `fig`

pp.save_plot(fig,  2018, "main_four", cfg2018; variant = "zero")


pp.save_plot(fig2, 2018, "no_forcing", cfg2018; variant = "zero")

pp.save_plot(fig3, 2018, "lula_bolsonaro", cfg2018; variant = "zero")




pp.@load "intermediate_data/group_metrics/group_metrics_2022.jld2" metrics


allgm_22 = copy(metrics )



pp.@load "intermediate_data/group_metrics/group_metrics_2006.jld2" metrics 


all_gm2006 = copy(metrics)

pp.@load "intermediate_data/group_metrics/group_metrics_2018.jld2" metrics


all_gm2018 = copy(metrics)


metrics = nothing 






scens = [(2006,"lula_alckmin"), (2018,"main_four"), (2022,"lula_bolsonaro")]

all_gm = Dict(2006=>all_gm2006, 2018=>all_gm2018, 2022=>allgm_22)


fig =pp.compare_demographic_across_scenarios(all_gm, f3, scens;
                                           demographic = "Ideology",
                                           variant     = :mice)


pp.save_plot(fig, 2022, "multiple_scenarios", f3[2022].cfg; variant = "mice")


scens = [(2006,"lula_alckmin"), (2022,"lula_bolsonaro")]



fig =pp.compare_demographic_across_scenarios(all_gm, f3, scens;
                                           demographic = "PT",
                                           variant     = :mice)


pp.save_plot(fig, 2022, "multiple_scenarios_PT", f3[2022].cfg; variant = "mice")




# Effective number of reversal preferences 


meas2022

meas2022["lula_bolsonaro"]

meas2022["lula_bolsonaro"][6]







tbl2022 = pp.enrp_table(Dict(2022 => meas2022), f3, 2022, "lula_bolsonaro";
                 variant = :mice, m_max = 7)





tbl2006 = pp.enrp_table(Dict(2006 => meas2006), f3, 2006, "lula_alckmin";
                 variant = :mice, m_max = 6)                 




tbl2018 = pp.enrp_table(Dict(2018 => meas2018), f3, 2018, "main_four";
                 variant = :mice, m_max = 7)        
                 
                 

tbl2018b = pp.enrp_table(Dict(2018 => meas2018), f3, 2018, "lula_bolsonaro";
                 variant = :mice, m_max = 7)              
                 
tbl2018c = pp.enrp_table(Dict(2018 => meas2018), f3, 2018, "no_forcing";
                 variant = :mice, m_max = 7)                    




tbl2018c = pp.hhi_table(Dict(2018 => meas2018), f3, 2018, "no_forcing";
                 variant = :mice, m_max = 7)    



tbl2006 = pp.polar_table(Dict(2006 => meas2006), f3, 2006, "lula_alckmin";
                       variant = :mice, m_max = 6)

tbl2018a = pp.polar_table(Dict(2018 => meas2018), f3, 2018, "main_four";
                        variant = :mice, m_max = 7)

tbl2018b = pp.polar_table(Dict(2018 => meas2018), f3, 2018, "lula_bolsonaro";
                        variant = :mice, m_max = 7)


tbl2018c = pp.polar_table(Dict(2018 => meas2018), f3, 2018, "no_forcing";
                        variant = :mice, m_max = 7)                        


tbl2022 = pp.polar_table(Dict(2022 => meas2022), f3, 2022, "lula_bolsonaro";
                       variant = :mice, m_max = 7)
