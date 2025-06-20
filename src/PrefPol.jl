module PrefPol

using Pkg

using CairoMakie
using CategoricalArrays
using Combinatorics
using DataFrames
using DataVoyager
using Dates 
import Impute 
using LaTeXStrings
using Pkg
using Random
using RCall                   # Enables calling R from Julia
using Statistics
using StatsBase
using CategoricalArrays
using PrettyTables
import  ProgressMeter as pm 
using Printf
using KernelDensity
using PyCall
import PyPlot as plt 
import Colors
using TextWrap
import OrderedCollections
using OrderedCollections: OrderedDict
using Logging

using JLD2, Arrow, TOML          # ← only for the helpers at bottom
# the measure with all the 13 candidates will probably be utterly meaningless!! 

project_root = dirname(Pkg.project().path)

dF_2022 = project_root * "/data/datafolha_vespera_2022_04780/04780/04780.SAV"
dF_2018 = project_root * "/data/04619/04619.SAV"




eseb_22 = project_root * "/data/04810/04810.sav"

eseb_18  = project_root *"/data/eseb_2018/04622/04622.sav"
eseb_06 = project_root * "/data/02489/1_02489.sav"


CANDIDATOS_eseb2022 = [
        "CIRO_GOMES", "BOLSONARO", "ALVARO_DIAS", "ARTHUR_LIRA", "LULA",
        "GERALDO_ALCKMIN", "GILBERTO_KASSAB", "EDUARDO_LEITE", "BOULOS",
        "MARINA_SILVA", "TARCISIO_DE_FREITAS", "LUCIANO_BIVAR", "SIMONE_TEBET"
    ]
include("preprocessing_general.jl")
include("preprocessing_specific.jl")
include("polarization_measures.jl")
include("mallows_play.jl")
include("plotting_bts.jl")
include("summary_measures.jl")
include("pipeline.jl")

#include("newplotting.jl")

export project_root, eseb_22, CANDIDATOS_eseb2022




end # module PrefPol
