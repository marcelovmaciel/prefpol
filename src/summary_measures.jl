
"""
    summarize_measures(measure_dict::Dict)

Given a nested dictionary of the form:
  :measure => :variant => Vector{Float64}
computes summary stats (mean, median, mode, min, max, Q1, Q3) for each measure/variant.

Returns a long-format DataFrame with columns:
  :measure, :variant, :stat, :value
"""
function summarize_measures(measure_dict::Dict)
    rows = NamedTuple[]

    for (measure, variant_dict) in measure_dict
        for (variant, vec) in variant_dict
            # Compute stats
            μ     = mean(vec)
            med   = median(vec)
            q1, q3 = quantile(vec, [0.25, 0.75])
            minv  = minimum(vec)
            maxv  = maximum(vec)
            modev = mode(vec)

            # Push all stats as separate rows
            append!(rows, [
                (measure = String(measure), variant = String(variant), stat = "mean",   value = μ),
                (measure = String(measure), variant = String(variant), stat = "median", value = med),


                (measure = String(measure), variant = String(variant), stat = "min",    value = minv),
                (measure = String(measure), variant = String(variant), stat = "max",    value = maxv),
            ])
        end
         #(measure = String(measure), variant = String(variant), stat = "mode",   value = modev),
          #      (measure = String(measure), variant = String(variant), stat = "q1",     value = q1),
           #     (measure = String(measure), variant = String(variant), stat = "q3",     value = q3),
               
    end

    return DataFrame(rows)
end



function summarize_measures_labeled(measure_dict::Dict, labels::Dict)
    rows = NamedTuple[]

    for (measure, variant_dict) in measure_dict
        # Lookup display label
        display_label = get(labels, measure, string(measure))

        for (variant, vec) in variant_dict
            # Compute stats
            μ     = round(mean(vec), digits = 3)
            med   = round(median(vec), digits = 3)
            q1, q3 = quantile(vec, [0.25, 0.75])
            minv  = round(minimum(vec), digits = 3)
            maxv  = round(maximum(vec), digits = 3)
            modev = mode(vec)

            # Push all stats as separate rows
            append!(rows, [
                (measure = display_label, variant = String(variant), stat = "mean",   value = μ),
                (measure = display_label, variant = String(variant), stat = "median", value = med),
            
                (measure = display_label, variant = String(variant), stat = "min",    value = minv),
                (measure = display_label, variant = String(variant), stat = "max",    value = maxv),
            ])
        end
    end

    return DataFrame(rows)
end




"""
    split_variant_tables(summary_df::DataFrame)

Given the long-format summary DataFrame (columns: :measure, :variant, :stat, :value),
returns a Dict with keys being variant names (Strings), and values DataFrames with:

| measure | mean | median | min | max |
"""
function split_variant_tables(summary_df::DataFrame)
    result = Dict{String, DataFrame}()

    variants = unique(summary_df.variant)

    for v in variants
        subdf = Base.filter(:variant => ==(v), summary_df)

        # Pivot stat column to wide format
        wide_df = unstack(subdf, :measure, :stat, :value)

        result[v] = wide_df
    end

    return result
end





"""
    save_variant_tables_latex(tables::Dict{String, DataFrame}, basename::String)

For each entry in `tables` (e.g., "zero", "random", "mice"), write a LaTeX file named
    "<basename>_<variant>.tex"
to the current directory.
"""
function save_variant_tables_latex(tables::Dict{String, DataFrame}, basename::String)
    for (variant, df) in tables
        filename = "$(basename)_$(variant).tex"
        open(filename, "w") do io
            pretty_table(io, df;
                backend = Val(:latex),
                tf = tf_latex_booktabs,   # ← NO COMMA HERE (this was the bug)
                alignment = :c,
                header_alignment = :c,
                title = "Imputation Variant: $(variant)")
        end
        println("Saved LaTeX table to: $filename")
    end
end




function save_sidebyside_latex(tables::Dict{String,DataFrame},
                               variants::Vector{String},
                               filename::String)

    stats     = [:mean, :median, :min, :max]               # assumed order
    stat_str  = ["mean", "median", "min", "max"]

    # 1 ── build the combined DataFrame
    base      = select(tables[variants[1]], :measure)      # first column

    for v in variants
        dfv = select(tables[v], Not(:measure))             # drop measure col
        # rename each stat column → "<variant>_<stat>"
        rename!(dfv, Symbol.(string.(v, "_", stat_str)))
        base = hcat(base, dfv)
    end

    # 2 ── build PrettyTables headers with spanners
    # first row: top spanners ("", 4 cols per variant)
    top_hdr = vcat([""; repeat(variants, inner=length(stats))]...)
    # second row: sub-headers for stats
    sub_hdr = vcat(["measure"; repeat(stat_str, outer=length(variants))]...)

    open(filename, "w") do io
        pretty_table(io, base;
            backend = Val(:latex),
            tf      = tf_latex_booktabs,    # booktabs styl
            header  = (top_hdr, sub_hdr),
            alignment = :c,
            header_alignment = :c,
            title   = "Summary Statistics by Imputation Variant")
    end
    println("Saved side-by-side LaTeX table → $filename")
end


