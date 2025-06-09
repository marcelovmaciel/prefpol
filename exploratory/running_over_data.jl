# Refactored ESEB Analysis Code
# =====================================

using Revise

using PrefPol

import PrefPol as pp

using RCall


using LaTeXStrings

using PyCall 
using DataFrames, Statistics

using Colors, ColorSchemes

# Configuration Module
module ESEBConfig
    # Election year configurations
    const ELECTIONS = Dict(
        2022 => Dict(
            :candidates => [
                "CIRO_GOMES", "BOLSONARO", "ALVARO_DIAS", "ARTHUR_LIRA", "LULA",
                "GERALDO_ALCKMIN", "GILBERTO_KASSAB", "EDUARDO_LEITE", "BOULOS",
                "MARINA_SILVA", "TARCISIO_DE_FREITAS", "LUCIANO_BIVAR", "SIMONE_TEBET"
            ],
            :demographics => ["Sex", "Religion", "Race", "Ideology", "PT", "Abortion"],
            :data_loader => :load_and_prepare_scores_df,
            :max_candidates => 7,
            :forced_scenarios => [
                Dict(:name => "lula_bolsonaro", :candidates => ["LULA", "BOLSONARO"])
            ]
        ),
        2018 => Dict(
            :candidates => [
                "Ciro_Gomes", "Manuela", "Guilherme_Boulos", "Marina_Silva",
                "Fernando_Haddad", "Henrique_Meirelles", "Jair_Bolsonaro",
                "Geraldo_Alckmin", "João_Amoêdo", "Lula", "Alvaro_Dias",
                "João_Goulart_Filho", "Cabo_Daciolo", "Rodrigo_Maia", "Eymael",
                "Vera", "Aécio_Neves", "Dilma_Rousseff", "Romero_Jucá",
                "Renan_Calheiros", "Michel_Temer"
            ],
            :demographics => ["Sex", "Religion", "Race", "Ideology"],
            :data_loader => :load_and_prepare_e2018,
            :max_candidates => 7,
            :forced_scenarios => [
                Dict(:name => "no_forcing", :candidates => String[]),
                Dict(:name => "lula_bolsonaro", :candidates => ["Lula", "Jair_Bolsonaro", "Ciro_Gomes", "Geraldo_Alckmin"]),
                Dict(:name => "main_four", :candidates => ["Fernando_Haddad", "Jair_Bolsonaro", "Ciro_Gomes", "Geraldo_Alckmin"])
            ]
        ),
        2006 => Dict(
            :candidates => [
                "Lula", "Geraldo_Alckmin", "Heloísa_Helena",
                "Cristóvam_Buarque", "Aécio_Neves", "José_Serra"
            ],
            :demographics => ["Sex", "Ideology", "PT"],
            :data_loader => :load_and_prepare_e2006,
            :max_candidates => 6,
            :forced_scenarios => [
                Dict(:name => "no_forcing", :candidates => String[]),
                Dict(:name => "lula_alckmin", :candidates => ["Lula", "Geraldo_Alckmin"])
            ]
        )
    )
    
    # Analysis parameters
    const DEFAULT_BOOTSTRAP_SAMPLES = Dict(
        :production => 1500,
        :development => 10,
        :testing => 5
    )
    
    const M_VALUES_RANGE = 2:7
end

# Analysis Results Structure
mutable struct ElectionAnalysisResults
    year::Int
    scenario_name::String
    candidate_set::Vector{String}
    measures_over_m::Any
    group_stats::Any
    global_plot::Any
    group_plot::Any
    
    # Constructor for partial initialization
    function ElectionAnalysisResults(year, scenario_name, candidate_set)
        new(year, scenario_name, candidate_set, nothing, nothing, nothing, nothing)
    end
    
    # Constructor for complete initialization
    function ElectionAnalysisResults(year, scenario_name, candidate_set, measures_over_m, group_stats, global_plot, group_plot)
        new(year, scenario_name, candidate_set, measures_over_m, group_stats, global_plot, group_plot)
    end
end

# Main Analysis Engine
struct ESEBAnalyzer
    config::Dict
    bootstrap_mode::Symbol  # :production, :development, :testing
    output_dir::String
    
    function ESEBAnalyzer(bootstrap_mode=:development, output_dir="exploratory/imgs")
        new(ESEBConfig.ELECTIONS, bootstrap_mode, output_dir)
    end
end

# Core analysis functions
function load_election_data(analyzer::ESEBAnalyzer, year::Int)
    election_config = analyzer.config[year]
    data_loader = election_config[:data_loader]
    candidates = election_config[:candidates]
    
    if data_loader == :load_and_prepare_scores_df
        return pp.load_and_prepare_scores_df(eseb_22)
    elseif data_loader == :load_and_prepare_e2018
        return pp.load_and_prepare_e2018(pp.eseb_18, CANDIDATOS=candidates)
    elseif data_loader == :load_and_prepare_e2006
        return pp.load_and_prepare_e2006(pp.eseb_06, candidates=candidates)
    else
        error("Unknown data loader: $data_loader")
    end
end

function get_bootstrap_samples(analyzer::ESEBAnalyzer)::Int
    return ESEBConfig.DEFAULT_BOOTSTRAP_SAMPLES[analyzer.bootstrap_mode]
end

function compute_candidate_sets(analyzer::ESEBAnalyzer, data, year::Int)
    election_config = analyzer.config[year]
    candidate_sets = Dict{String, Vector{String}}()
    
    for scenario in election_config[:forced_scenarios]
        scenario_name = scenario[:name]
        forced_candidates = scenario[:candidates]
        
        if isempty(forced_candidates)
            candidate_set = pp.compute_candidate_set(
                data;
                candidate_cols = election_config[:candidates],
                m = election_config[:max_candidates]
            )
        else
            candidate_set = pp.compute_candidate_set(
                data;
                candidate_cols = election_config[:candidates],
                m = election_config[:max_candidates],
                force_include = forced_candidates
            )
        end
        
        candidate_sets[scenario_name] = candidate_set
    end
    
    return candidate_sets
end

function analyze_global_measures(analyzer::ESEBAnalyzer, data, year::Int, forced_candidates::Vector{String}=String[])
    election_config = analyzer.config[year]
    m_values = collect(ESEBConfig.M_VALUES_RANGE)[1:election_config[:max_candidates]-1]
    bootstrap_samples = get_bootstrap_samples(analyzer)
    
    kwargs = Dict(
        :n_bootstrap => bootstrap_samples,
        :candidate_cols => election_config[:candidates],
        :demographics => election_config[:demographics]
    )
    
    if !isempty(forced_candidates)
        kwargs[:force_include] = forced_candidates
    end
    
    return pp.compute_measures_over_alternatives(data, m_values; kwargs...)
end

function analyze_group_measures(analyzer::ESEBAnalyzer, data, year::Int, forced_candidates::Vector{String}=String[])
    election_config = analyzer.config[year]
    m_values = collect(ESEBConfig.M_VALUES_RANGE)[1:election_config[:max_candidates]-1]
    bootstrap_samples = get_bootstrap_samples(analyzer)
    
    kwargs = Dict(
        :candidate_cols => election_config[:candidates],
        :demographics => election_config[:demographics],
        :m_values => m_values,
        :n_bootstrap => bootstrap_samples,
        :force_include => forced_candidates
    )
    
    return pp.compute_stats_over_m(data; kwargs...)
end

function create_plots(analyzer::ESEBAnalyzer, results::ElectionAnalysisResults, measures_data, group_stats; variants = ["zero", "random", "mice"])
    election_config = analyzer.config[results.year]
    candidate_label = pp.describe_candidate_set(results.candidate_set)
    
    # Create global measures plot
    global_plot = pp.lines_alt_by_variant(
        measures_data, 
        candidate_label = candidate_label, 
        year = results.year,
        variants = variants
    )

    
    # Create group measures plot
    maxcols = results.year == 2006 ? 1 : 2
    group_plot = pp.lines_group_measures_over_m(
        group_stats,
        demographics = Symbol.(election_config[:demographics]),
        candidate_label = candidate_label,
        year = results.year,
        maxcols = maxcols,
        variants = variants
    )
    
    return global_plot, group_plot
end

function save_plots(analyzer::ESEBAnalyzer, results::ElectionAnalysisResults, global_plot, group_plot)
    bootstrap_suffix = "$(get_bootstrap_samples(analyzer))bts"
    
    global_filename = "alt_by_variant_lines$(results.year)_$(bootstrap_suffix)_$(results.scenario_name).png"
    group_filename = "cgd_lines$(results.year)_$(bootstrap_suffix)_$(results.scenario_name).png"
    
    global_path = joinpath(analyzer.output_dir, global_filename)
    group_path = joinpath(analyzer.output_dir, group_filename)
    
    pp.save(global_path, global_plot; px_per_unit=3)
    pp.save(group_path, group_plot; px_per_unit=3)
    
    @info "Saved plots for $(results.year) $(results.scenario_name): $global_filename, $group_filename"
end

function analyze_single_scenario(analyzer::ESEBAnalyzer, data, year::Int, scenario_name::String, candidate_set::Vector{String}, variants_to_plot = ["mice"])
    election_config = analyzer.config[year]
    
    # Find the forced candidates for this scenario
    scenario_config = findfirst(s -> s[:name] == scenario_name, election_config[:forced_scenarios])
    forced_candidates = scenario_config !== nothing ? election_config[:forced_scenarios][scenario_config][:candidates] : String[]
    
    @info "Analyzing $year - $scenario_name scenario..."
    
    # Initialize results
    results = ElectionAnalysisResults(year, scenario_name, candidate_set)
    
    # Analyze global measures
    @info "  Computing global measures..."
    measures_data = analyze_global_measures(analyzer, data, year, forced_candidates)
    
    # Analyze group measures  
    @info "  Computing group measures..."
    group_stats = analyze_group_measures(analyzer, data, year, forced_candidates)
    
    # Create plots
    @info "  Creating plots..."
    global_plot, group_plot = create_plots(analyzer, results, measures_data, group_stats, variants = variants_to_plot)
    
    # Save plots
    save_plots(analyzer, results, global_plot, group_plot)
    
    # Update the results with computed data
    results.measures_over_m = measures_data
    results.group_stats = group_stats
    results.global_plot = global_plot
    results.group_plot = group_plot
    
    return results
end

function analyze_election_year(analyzer::ESEBAnalyzer, year::Int)
    @info "Starting analysis for election year $year"
    
    # Load data
    @info "Loading data for $year..."
    data = load_election_data(analyzer, year)
    
    # Compute candidate sets for all scenarios
    @info "Computing candidate sets..."
    candidate_sets = compute_candidate_sets(analyzer, data, year)
    
    # Analyze each scenario
    results = Dict{String, ElectionAnalysisResults}()
    for (scenario_name, candidate_set) in candidate_sets
        results[scenario_name] = analyze_single_scenario(analyzer, data, year, scenario_name, candidate_set)
    end
    
    @info "Completed analysis for $year"
    return results
end

function run_full_analysis(analyzer::ESEBAnalyzer, years::Vector{Int}=[2006, 2018, 2022])
    @info "Starting full ESEB analysis for years: $years"
    @info "Bootstrap mode: $(analyzer.bootstrap_mode) ($(get_bootstrap_samples(analyzer)) samples)"
    
    all_results = Dict{Int, Dict{String, ElectionAnalysisResults}}()
    
    for year in years
        try
            all_results[year] = analyze_election_year(analyzer, year)
        catch e
            @error "Failed to analyze year $year: $e"
            rethrow(e)
        end
    end
    
    @info "Full analysis completed successfully!"
    return all_results
end

# Convenience functions for different use cases
function run_development_analysis(years::Vector{Int}=[2022])
    analyzer = ESEBAnalyzer(:development)
    return run_full_analysis(analyzer, years)
end

function run_production_analysis(years::Vector{Int}=[2006, 2018, 2022])
    analyzer = ESEBAnalyzer(:production)
    return run_full_analysis(analyzer, years)
end

function run_quick_test(year::Int=2022)
    analyzer = ESEBAnalyzer(:testing)
    return run_full_analysis(analyzer, [year])
end

# Main execution - replace the original messy code with this:
function main()
    # For development/testing - fast execution
    # results = run_development_analysis([2022])
    
    # For production - full analysis with high bootstrap samples
    results = run_production_analysis()
    
    return results
end

# Execute the analysis



#prod = run_production_analysis([2006, 2018, 2022])


# Save
#@save "prod15kbts.jld2" prod

# Load (as variable)
#@load "filename.jld2" my_object



fig_global_2018 = pp.compare_global_measures_v3(prod[2018];
                                            variants = ["mice"])


pp.save("exploratory/imgs/globalvars_2018Scenarios.png", fig_global_2018,  px_per_unit=3)


fig_global_2006 = pp.compare_global_measures_v3(prod[2006];
                                            variants = ["mice"])

pp.save("exploratory/imgs/globalvars_2006Scenarios.png", fig_global_2006,  px_per_unit=3)


#= fig = pp.compare_demographic_across_scenarios(
          all_results,
          [(2006, "lula_alckmin"),
           (2018, "main_four"),
           (2022, "lula_bolsonaro")];
          demographic = "Ideology")  =#


figSex = pp.compare_demographic_across_scenariosy(
          prod,
          [(2006, "lula_alckmin"),
           (2018, "main_four"),
           (2022, "lula_bolsonaro")];
          demographic = "Sex")           

pp.save("exploratory/imgs/sex_comparison.png", figSex; px_per_unit=3)          




figPT = pp.compare_demographic_across_scenariosy(
          prod,
          [(2006, "lula_alckmin"),
           (2022, "lula_bolsonaro")];
          demographic = "PT")    

pp.save("exploratory/imgs/pt_comparison.png", figPT; px_per_unit=3)



fig_race = pp.compare_demographic_across_scenariosy(
          prod,
          [(2018, "main_four"),
           (2018, "lula_bolsonaro"),
           (2022, "lula_bolsonaro")];
          demographic = "Race")  
          
pp.save("exploratory/imgs/race_comparison.png", fig_race; px_per_unit=3)  



fig_religion = pp.compare_demographic_across_scenariosy(
          prod,
          [(2018, "main_four"),
           (2018, "lula_bolsonaro"),
           (2022, "lula_bolsonaro")];
          demographic = "Religion")  
          
pp.save("exploratory/imgs/religion_comparison.png", fig_religion; px_per_unit=3) 



fig_ideology = pp.compare_demographic_across_scenariosy(
          prod,
          [(2006, "lula_alckmin"),
           (2018, "main_four"),
           (2022, "lula_bolsonaro")];
          demographic = "Ideology")  
          
pp.save("exploratory/imgs/ideology_comparison.png", fig_ideology; px_per_unit=3) 



keys(prod[2022]["lula_bolsonaro"])

prod[2022]["lula_bolsonaro"]