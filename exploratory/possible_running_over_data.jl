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
    
    const M_VALUES_RANGE = 2:7
end

# Bootstrap Configuration Struct
struct BootstrapConfig
    year::Int
    scenario_name::String
    candidate_set::Vector{String}
    forced_candidates::Vector{String}
    n_bootstrap::Int
    m_values::Vector{Int}
    demographics::Vector{String}
    candidate_cols::Vector{String}
    
    function BootstrapConfig(
        year::Int, 
        scenario_name::String, 
        candidate_set::Vector{String},
        forced_candidates::Vector{String} = String[];
        n_bootstrap::Int = 1500,
        m_values::Vector{Int} = collect(2:6),
        demographics::Vector{String} = String[],
        candidate_cols::Vector{String} = String[]
    )
        # Get election config for validation and defaults
        election_config = ESEBConfig.ELECTIONS[year]
        
        # Use defaults from election config if not provided
        final_demographics = isempty(demographics) ? election_config[:demographics] : demographics
        final_candidate_cols = isempty(candidate_cols) ? election_config[:candidates] : candidate_cols
        final_m_values = isempty(m_values) ? collect(2:election_config[:max_candidates]) : m_values
        
        new(year, scenario_name, candidate_set, forced_candidates, n_bootstrap, 
            final_m_values, final_demographics, final_candidate_cols)
    end
end



struct AnalysisResults
    config::BootstrapConfig
    global_measures::Any
    group_measures::Any
    bootstraps::Dict{Symbol,Vector{DataFrame}}      # NEW field

    function AnalysisResults(config::BootstrapConfig,
                             global_measures,
                             group_measures,
                             bootstraps::Dict{Symbol,Vector{DataFrame}})
        new(config, global_measures, group_measures, bootstraps)
    end
end
# Plot Results Structure
struct PlotResults
    analysis_results::AnalysisResults
    global_plot::Any
    group_plot::Any
    
    function PlotResults(analysis_results::AnalysisResults, global_plot, group_plot)
        new(analysis_results, global_plot, group_plot)
    end
end

# Data Loading Module
struct DataLoader
    config::Dict
    
    function DataLoader()
        new(ESEBConfig.ELECTIONS)
    end
end

function load_election_data(loader::DataLoader, year::Int)
    election_config = loader.config[year]
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

function compute_candidate_sets(loader::DataLoader, data, year::Int)
    election_config = loader.config[year]
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

# Analysis Engine - Only computes measures
struct MeasureCalculator
    data_loader::DataLoader
    
    function MeasureCalculator()
        new(DataLoader())
    end
end

function compute_global_measures(calculator::MeasureCalculator, data, config::BootstrapConfig)
    @info "Computing global measures for $(config.year) $(config.scenario_name)..."
    
    kwargs = Dict(
        :n_bootstrap => config.n_bootstrap,
        :candidate_cols => config.candidate_cols,
        :demographics => config.demographics
    )
    
    if !isempty(config.forced_candidates)
        kwargs[:force_include] = config.forced_candidates
    end
    
    return pp.compute_measures_over_alternatives(data, config.m_values; kwargs...)
end

function compute_group_measures(calculator::MeasureCalculator, data, config::BootstrapConfig)
    @info "Computing group measures for $(config.year) $(config.scenario_name)..."
    
    kwargs = Dict(
        :candidate_cols => config.candidate_cols,
        :demographics => config.demographics,
        :m_values => config.m_values,
        :n_bootstrap => config.n_bootstrap,
        :force_include => config.forced_candidates
    )
    
    return pp.compute_stats_over_m(data; kwargs...)
end



"""
    generate_bootstraps(df, cfg) → Dict{Symbol,Vector{DataFrame}}

• Applies the three imputation variants (zero, random, mice)  
• Draws `cfg.n_bootstrap` weighted replicates for each variant  
• Returns the `Dict(:zero=>[..], :random=>[..], :mice=>[..])`
"""
function generate_bootstraps(df::DataFrame, cfg::BootstrapConfig)
    imps = pp.imputation_variants(df,               # ← the imputation you already wrote
                               cfg.candidate_cols,
                               cfg.demographics;
                               most_known_candidates = cfg.candidate_set)

    # survey weight column is usually called :peso or :weight – adapt if needed
    w = hasproperty(df, :peso) ? Vector{Float64}(df.peso) :
        hasproperty(df, :weight) ? Vector{Float64}(df.weight) :
        ones(Float64, nrow(df))            # fallback: equal weights

    return pp.bootstrap_variants(imps; weights = w, B = cfg.n_bootstrap)
end


function analyze_scenario(calculator::MeasureCalculator,
                          data::DataFrame,
                          cfg::BootstrapConfig)

    @info "Analyzing scenario: $(cfg.year) – $(cfg.scenario_name)"

    # ─── ❶ generate & store replicates ─────────────────────────
    boots = generate_bootstraps(data, cfg)

    # ─── ❷ your existing point-estimates ───────────────────────
    global_meas = compute_global_measures(calculator, data, cfg)
    group_meas  = compute_group_measures(calculator,  data, cfg)

    # ─── ❸ bundle everything ──────────────────────────────────
    return AnalysisResults(cfg, global_meas, group_meas, boots)
end



# Plotting Engine - Separate from calculation
struct PlotGenerator
    output_dir::String
    
    function PlotGenerator(output_dir::String = "exploratory/imgs")
        new(output_dir)
    end
end

function create_global_plot(generator::PlotGenerator, results::AnalysisResults; variants = ["zero", "random", "mice"])
    config = results.config
    candidate_label = pp.describe_candidate_set(config.candidate_set)
    
    return pp.lines_alt_by_variant(
        results.global_measures, 
        candidate_label = candidate_label, 
        year = config.year,
        variants = variants
    )
end

function create_group_plot(generator::PlotGenerator, results::AnalysisResults; variants = ["zero", "random", "mice"])
    config = results.config
    candidate_label = pp.describe_candidate_set(config.candidate_set)
    maxcols = config.year == 2006 ? 1 : 2
    
    return pp.lines_group_measures_over_m(
        results.group_measures,
        demographics = Symbol.(config.demographics),
        candidate_label = candidate_label,
        year = config.year,
        maxcols = maxcols,
        variants = variants
    )
end

function create_plots(generator::PlotGenerator, results::AnalysisResults; variants = ["zero", "random", "mice"])
    @info "Creating plots for $(results.config.year) $(results.config.scenario_name)..."
    
    global_plot = create_global_plot(generator, results, variants = variants)
    group_plot = create_group_plot(generator, results, variants = variants)
    
    return PlotResults(results, global_plot, group_plot)
end

function save_plots(generator::PlotGenerator, plot_results::PlotResults)
    config = plot_results.analysis_results.config
    bootstrap_suffix = "$(config.n_bootstrap)bts"
    
    global_filename = "alt_by_variant_lines$(config.year)_$(bootstrap_suffix)_$(config.scenario_name).png"
    group_filename = "cgd_lines$(config.year)_$(bootstrap_suffix)_$(config.scenario_name).png"
    
    global_path = joinpath(generator.output_dir, global_filename)
    group_path = joinpath(generator.output_dir, group_filename)
    
    pp.save(global_path, plot_results.global_plot; px_per_unit=3)
    pp.save(group_path, plot_results.group_plot; px_per_unit=3)
    
    @info "Saved plots: $global_filename, $group_filename"
end

# Configuration Builder Functions
function create_bootstrap_configs(year::Int, n_bootstrap::Int = 1500)
    loader = DataLoader()
    data = load_election_data(loader, year)
    candidate_sets = compute_candidate_sets(loader, data, year)
    
    configs = BootstrapConfig[]
    election_config = ESEBConfig.ELECTIONS[year]
    
    for scenario in election_config[:forced_scenarios]
        scenario_name = scenario[:name]
        candidate_set = candidate_sets[scenario_name]
        forced_candidates = scenario[:candidates]
        
        config = BootstrapConfig(
            year, 
            scenario_name, 
            candidate_set, 
            forced_candidates;
            n_bootstrap = n_bootstrap
        )
        
        push!(configs, config)
    end
    
    return configs, data
end

# High-level Workflow Functions
function run_analysis_pipeline(configs::Vector{BootstrapConfig}, data; 
                             create_plots_flag::Bool = true, 
                             save_plots_flag::Bool = true,
                             plot_variants = ["mice"])
    calculator = MeasureCalculator()
    plot_generator = PlotGenerator()
    
    analysis_results = AnalysisResults[]
    plot_results = PlotResults[]
    
    # Calculate measures for all configurations
    @info "Computing measures for $(length(configs)) configurations..."
    for config in configs
        result = analyze_scenario(calculator, data, config)
        push!(analysis_results, result)
    end
    
    # Create plots if requested
    if create_plots_flag
        @info "Creating plots..."
        for result in analysis_results
            plot_result = create_plots(plot_generator, result, variants = plot_variants)
            push!(plot_results, plot_result)
            
            if save_plots_flag
                save_plots(plot_generator, plot_result)
            end
        end
    end
    
    return analysis_results, plot_results
end

function analyze_year(year::Int; n_bootstrap::Int = 1500, create_plots_flag::Bool = true, save_plots_flag::Bool = true)
    @info "Starting analysis for year $year with $n_bootstrap bootstrap samples"
    
    configs, data = create_bootstrap_configs(year, n_bootstrap)
    return run_analysis_pipeline(configs, data; 
                               create_plots_flag = create_plots_flag, 
                               save_plots_flag = save_plots_flag)
end

function analyze_multiple_years(years::Vector{Int}; n_bootstrap::Int = 1500, create_plots_flag::Bool = true, save_plots_flag::Bool = true)
    all_analysis_results = Dict{Int, Vector{AnalysisResults}}()
    all_plot_results = Dict{Int, Vector{PlotResults}}()
    
    for year in years
        analysis_results, plot_results = analyze_year(year; 
                                                    n_bootstrap = n_bootstrap,
                                                    create_plots_flag = create_plots_flag,
                                                    save_plots_flag = save_plots_flag)
        all_analysis_results[year] = analysis_results
        all_plot_results[year] = plot_results
    end
    
    return all_analysis_results, all_plot_results
end

# Convenience functions for different use cases
function run_development_analysis(years::Vector{Int} = [2022])
    return analyze_multiple_years(years; n_bootstrap = 10, create_plots_flag = true, save_plots_flag = false)
end

function run_production_analysis(years::Vector{Int} = [2006, 2018, 2022])
    return analyze_multiple_years(years; n_bootstrap = 1500, create_plots_flag = true, save_plots_flag = true)
end

function run_quick_test(year::Int = 2022)
    return analyze_year(year; n_bootstrap = 5, create_plots_flag = false, save_plots_flag = false)
end

# Analysis-only functions (no plotting)
function compute_measures_only(years::Vector{Int}; n_bootstrap::Int = 1500)
    return analyze_multiple_years(years; n_bootstrap = n_bootstrap, create_plots_flag = false, save_plots_flag = false)
end

# Plotting-only functions (from existing results)
function create_plots_from_results(analysis_results::Vector{AnalysisResults}; 
                                 variants = ["mice"], 
                                 save_plots_flag::Bool = true)
    plot_generator = PlotGenerator()
    plot_results = PlotResults[]
    
    for result in analysis_results
        plot_result = create_plots(plot_generator, result, variants = variants)
        push!(plot_results, plot_result)
        
        if save_plots_flag
            save_plots(plot_generator, plot_result)
        end
    end
    
    return plot_results
end

# Main execution function
function main()
    # Example usage:
    
    # 1. For development - fast execution, no saving
    # analysis_results, plot_results = run_development_analysis([2022])
    
    # 2. For production - full analysis with plots
    # analysis_results, plot_results = run_production_analysis()
    
    # 3. Compute measures only (for later plotting)
    # analysis_results, _ = compute_measures_only([2022])
    
    # 4. Create plots from existing analysis results
    # plot_results = create_plots_from_results(analysis_results[2022])
    
    # Default: production analysis
    return run_production_analysis()
end


configs, data = create_bootstrap_configs(2022, 1500)

function compute_bootstraps_only(year::Int; n_bootstrap::Int = 1500)
    # ❶ Data + configs
    configs, data = create_bootstrap_configs(year, n_bootstrap)
    results = Dict{String,Dict}()

    # ❷ For each scenario: generate replicates, optionally save them
    for cfg in configs
        br = generate_bootstraps(data, cfg)                    # ← fast
        results[cfg.scenario_name] = br
        # uncomment next line if you also want an on-disk copy
        # save_bootstraps("exploratory/imgs/bootstraps_$(year)_$(cfg.scenario_name).jld2", br)
    end
    return results
end


boots_main = compute_bootstraps_only(2018; n_bootstrap = 10)

boots_main["main_four"][:mice]
