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



elec06 = pp.load_election_cfg("config/2006.toml") 

elec18 = pp.load_election_cfg("config/2018.toml") 

elec22 = pp.load_election_cfg("config/2022.toml") 



df_raw06 = pp.load_spss_file(elec06.data_file)

df_raw06.eseb16a |> pp.proportionmap 


letters = ['a','b','c','d','e','f']

 function build_letter_column_symbols(base::AbstractString, letters::Vector{Char})
        return [Symbol(base * string(c)) for c in letters]
end

rename!(df_raw06, Dict(zip(build_letter_column_symbols("eseb16", letters), elec06.candidates)))

df_raw06[!, elec06.candidates]


pairs = (11.0 => 99.0, 77.0 => 99.0)
for col in elec06.candidates
    replace!(df_raw06[!, col], pairs...)
end



pp.pp_proportions(df_raw06, elec06.candidates)


df_raw06.peso = df_raw06.peso_1

df_raw06.Sex = pp.categorical(df_raw06.SEXO)

df_raw06.eseb15a |> pp.proportionmap

pp.pp_proportions(df_raw06, [:eseb15a])


replace!(df_raw06.eseb15a, pairs...)

pp.pp_proportions(df_raw06, [:eseb15a])



replace!(x -> ismissing(x) ? x : x < 5 ? 0.0 : x <= 10 ? 1.0 : 99.0,
             df_raw06[!, :eseb15a])

pp.pp_proportions(df_raw06, [:eseb15a])



pp.pp_proportions(df_raw06, [:eseb19])



pairs = (66.0 => 99.0, 77.0 => 99.0)

replace!(df_raw06[!, :eseb19], pairs...)

pp.pp_proportions(df_raw06, [:eseb19])


function recode19(x)
    if ismissing(x)
        return x                     # keep missing
    elseif x <= 3                    # 0–3  → –1
        return -1
    elseif x <= 6                    # 4–6  → 0
        return 0
    elseif x <= 10                   # 7–10 → 1
        return 1
    else
        return x                     # 99 (or anything >10) unchanged
    end
end

df_raw06.eseb19 .= recode19.(df_raw06.eseb19)   # broadcast + in-place assign


pp.pp_proportions(df_raw06, [:eseb19])





df_e18 = pp.load_spss_file(elec18.data_file)



function build_column_symbols(base::AbstractString, n::Integer;
                                minwidth::Int = 2)
        width = max(minwidth, length(string(n)))          # e.g. n = 21  → width = 2
        return [Symbol(base * @sprintf("%0*d", width, i)) for i in 1:n]
end


elec18.candidates

rename!(df_e18, Dict(zip(build_column_symbols("Q16", 21), elec18.candidates)))

pp.pp_proportions(df_e18, elec18.candidates)

pairs = (96.0 => 99.0, 97.0 => 99.0, 98.0 => 99.0)

for col in elec18.candidates
        replace!(df_e18[!, col], pairs...)
end


pp.pp_proportions(df_e18, elec18.candidates)

pp.pp_proportions(df_e18, [:D10])


replace!(x -> x in (97.) ? 96.0 : x, df_e18.D10)

replace!(x -> x in (98.) ? 99.0 : x, df_e18.D10)

df_e18.D10 = pp.categorical(df_e18.D10)
df_e18.Religion = df_e18.D10



pp.pp_proportions(df_e18, [:D10])

df_e18.Sex = pp.categorical(df_e18.D2_SEXO)


pp.pp_proportions(df_e18, [:Sex])

replace!(x -> x in (8., 9.) ? 9.0 : x, df_e18.D12A)

df_e18.Race = df_e18.D12A


pp.pp_proportions(df_e18, [:Race])

pp.pp_proportions(df_e18, [:Q18])



replace!(x -> x in (95.0, 97.0, 98.0) ? 99.0 : x, df_e18.Q18)

pp.pp_proportions(df_e18, [:Q18])



function recodeQ18(x)
        if ismissing(x)
            return x                     # keep missing
        elseif x <= 3                    # 0–3  → –1
            return -1
        elseif x <= 6                    # 4–6  → 0
            return 0
        elseif x <= 10                   # 7–10 → 1
            return 1
        else
            return x                     # 99 (or anything >10) unchanged
        end
end

df_e18.Ideology .= recodeQ18.(df_e18.Q18)   # broadcast + in-place assign



df_e18.Ideology= pp.categorical(df_e18.Ideology;
                ordered = true,
                levels  = [-1, 0, 1, 99])


pp.pp_proportions(df_e18, [:Ideology])

pp.pp_proportions(df_e18, [:Q1513])

pairs = (96.0 => 99.0, 97.0 => 99.0, 98.0 => 99.0)


replace!(df_e18.Q1513, pairs...)

pp.pp_proportions(df_e18, [:Q1513])


replace!(x -> ismissing(x) ? x : x < 5 ? 0.0 : x <= 10 ? 1.0 : 99.0,
             df_e18[!, :Q1513])

df_e18.PT = df_e18.Q1513



pp.pp_proportions(df_e18, [:PT])             


df_e22 = pp.load_spss_file(elec22.data_file)


build_column_symbols(base::String, n::Int) = [Symbol(base * string(i)) for i in 1:n]

rename!(df_e22, Dict(zip(build_column_symbols("Q17_", 13), elec22.candidates)))


pp.pp_proportions(df_e22, elec22.candidates)


pairs = (96.0 => 99.0, 97.0 => 99.0, 98.0 => 99.0)

for col in elec22.candidates
        replace!(df_e22[!, col], pairs...)
end


pp.pp_proportions(df_e22, elec22.candidates)


pp.pp_proportions(df_e22, [:D10])


# I recoded 99, 100, 101, 102 -> 95; 96, 97-> 96; and 97, 98 -> 99
replace!(x -> x in (99., 100., 101., 102.) ? 95.0 : x, df_e22.D10)
replace!(x -> x in (96., 97.) ? 96.0 : x, df_e22.D10)
replace!(x -> x in (97., 98.) ? 99.0 : x, df_e22.D10)

df_e22.Religion = pp.categorical(df_e22.D10)

pp.pp_proportions(df_e22, [:Religion])


    
# Recode D02 (sex)
df_e22.Sex = pp.categorical(df_e22.D02)

# Recode D12a (race)
replace!(x -> x in (97.0, 98.0) ? 99.0 : x, df_e22.D12a)
    
df_e22.Race = pp.categorical(df_e22.D12a)

pp.pp_proportions(df_e22, [:Race])


pp.pp_proportions(df_e22, [:Q19])

# 95,96,98 -> 99.


replace!(x -> x in (95.0, 96.0, 98.0) ? 99.0 : x, df_e22.Q19)

pp.pp_proportions(df_e22, [:Q19])


function recode19(x)
        if ismissing(x)
            return x                     # keep missing
        elseif x <= 3                    # 0–3  → –1
            return -1
        elseif x <= 6                    # 4–6  → 0
            return 0
        elseif x <= 10                   # 7–10 → 1
            return 1
        else
            return x                     # 99 (or anything >10) unchanged
        end
    end

df_e22.Q19 .= recode19.(df_e22.Q19)   # broadcast + in-place assign

df_e22.Ideology = pp.categorical(df_e22.Q19;
                ordered = true,
                levels  = [-1, 0, 1, 99])


pp.pp_proportions(df_e22, [:Ideology])    

pp.pp_proportions(df_e22, [:Q18_5])

# I recoded 95, 96,97,98 -> 99
replace!(x -> x in (95.0, 96.0, 97.0, 98.0) ? 99.0 : x, df_e22.Q18_5)


pp.pp_proportions(df_e22, [:Q18_5])

replace!(x -> ismissing(x) ? x : x < 5 ? 0.0 : x <= 10 ? 1.0 : 99.0,
                 df_e22[!, :Q18_5])

df_e22.PT = df_e22.Q18_5

pp.pp_proportions(df_e22, [:PT])

pp.pp_proportions(df_e22, [:Q31_7])

# 97,98 -> 99

replace!( x -> x in (97.0, 98.0) ? 99.0 : x, df_e22.Q31_7)


df_e22.Abortion = df_e22.Q31_7


pp.pp_proportions(df_e22, [:Abortion])
