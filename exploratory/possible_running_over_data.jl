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
elec = pp.load_election_cfg("config/2018.toml")   # full file into memory


# ===============================================
# from here I need to change 
bc   = pp.make_bootstrap_cfg(elec, "lula_bolsonaro"; n_bootstrap = 2000)
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