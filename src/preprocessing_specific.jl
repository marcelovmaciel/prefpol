function pp_proportions(df::DataFrame, cols)
    for col in cols
        v  = df[!, col]                  # keeps missings
        pm = proportionmap(v)            # Dict(value => share)
        N  = length(v)                   # total observations (incl. missings)

        println("\n" * "─"^40)
        @printf("%-15s │ %8s │ %s\n", string(col), "prop.", "count")
        println("─"^40)

        for (val, p) in sort(collect(pm); by = first)  # deterministic order
            @printf("%-15s │ %6.2f%% │ %d\n",
                    val, p * 100, Int(round(p * N)))
        end
    end
    println("─"^40)
end








function load_and_prepare_scores_df(data_path::String; candidates = CANDIDATOS_eseb2022)
    # Load SPSS file
    df_e22 = load_spss_file(data_path)

    # Metadata
    #PARTIDOS = ["PDT", "PL", "PODEMOS", "PP", "PT", "PSB", "PSD", "PSDB", "PSOL", "REDE", "REP", "UB", "MDB"]

   
    build_column_symbols(base::String, n::Int) = [Symbol(base * string(i)) for i in 1:n]

    rename!(df_e22, Dict(zip(build_column_symbols("Q17_", 13), candidates)))


    pairs = (96.0 => 99.0, 97.0 => 99.0, 98.0 => 99.0)

    for col in candidates
            replace!(df_e22[!, col], pairs...)
    end


    # Recode D10 (religion)
    replace!(x -> x in (99., 100., 101., 102.) ? 95.0 : x, df_e22.D10)
    replace!(x -> x in (96., 97.) ? 96.0 : x, df_e22.D10)
    replace!(x -> x in (97., 98.) ? 99.0 : x, df_e22.D10)

    df_e22.Religion = categorical(df_e22.D10)

    
    # Recode D02 (sex)
    df_e22.Sex = categorical(df_e22.D02)

    # Recode D12a (race)
    replace!(x -> x in (97.0, 98.0) ? 99.0 : x, df_e22.D12a)
    
    df_e22.Race = categorical(df_e22.D12a)
    

    # Q19 / Ideology 
    replace!(x -> x in (95.0, 96.0, 98.0) ? 99.0 : x, df_e22.Q19)
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
    df_e22.Ideology = categorical(df_e22.Q19;
                    ordered = true,
                    levels  = [-1, 0, 1, 99])
    
    # PT 
    replace!(x -> x in (95.0, 96.0, 97.0, 98.0) ? 99.0 : x, df_e22.Q18_5)
    replace!(x -> ismissing(x) ? x : x < 5 ? 0.0 : x <= 10 ? 1.0 : 99.0,
                 df_e22[!, :Q18_5])
    df_e22.PT = df_e22.Q18_5

   replace!( x -> x in (97.0, 98.0) ? 99.0 : x, df_e22.Q31_7)
   df_e22.Abortion = df_e22.Q31_7
    return df_e22
end



function load_and_prepare_e2006(df_path; candidates = candidates2006)
    df_e06 = load_spss_file(df_path)
    letters = ['a','b','c','d','e','f']

    function build_letter_column_symbols(base::AbstractString, letters::Vector{Char})
        return [Symbol(base * string(c)) for c in letters]
    end

    rename!(df_e06, Dict(zip(build_letter_column_symbols("eseb16", letters), candidates)))

    pairs = (11.0 => 99.0, 77.0 => 99.0)
    for col in candidates
        replace!(df_e06[!, col], pairs...)
    end


    df_e06.peso = df_e06.peso_1

    df_e06.Sex = categorical(df_e06.SEXO)

    replace!(df_e06.eseb15a, pairs...)
    replace!(x -> ismissing(x) ? x : x < 5 ? 0.0 : x <= 10 ? 1.0 : 99.0,
             df_e06[!, :eseb15a])

    df_e06.PT = df_e06.eseb15a
    
    pairs = (66.0 => 99.0, 77.0 => 99.0)

    replace!(df_e06[!, :eseb19], pairs...)
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
    df_e06.eseb19 .= recode19.(df_e06.eseb19)   # broadcast + in-place assign



    df_e06.Ideology= categorical(df_e06.eseb19;
                ordered = true,
                levels  = [-1, 0, 1, 99])

    return(df_e06)
end



function load_and_prepare_e2018(df_path; candidates = candidates2018)
   

    df_e18 = load_spss_file(df_path)

    function build_column_symbols(base::AbstractString, n::Integer;
                                minwidth::Int = 2)
        width = max(minwidth, length(string(n)))          # e.g. n = 21  → width = 2
        return [Symbol(base * @sprintf("%0*d", width, i)) for i in 1:n]
    end



    rename!(df_e18, Dict(zip(build_column_symbols("Q16", 21), candidates)))
    pairs = (96.0 => 99.0, 97.0 => 99.0, 98.0 => 99.0)
    for col in candidates
        replace!(df_e18[!, col], pairs...)
    end

    ##  religion =====================================================================================================================
    replace!(x -> x in (97.) ? 96.0 : x, df_e18.D10)
    replace!(x -> x in (98.) ? 99.0 : x, df_e18.D10)
    df_e18.D10 = categorical(df_e18.D10)
    df_e18.Religion = df_e18.D10

    # sex =====================================================================================================================
    df_e18.Sex = categorical(df_e18.D2_SEXO)

    # Race ===================================================================================================================== 
    replace!(x -> x in (8., 9.) ? 9.0 : x, df_e18.D12A)

    df_e18.Race = df_e18.D12A


    # ideology ==============================================================================
    replace!(x -> x in (95.0, 97.0, 98.0) ? 99.0 : x, df_e18.Q18)




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



    df_e18.Ideology= categorical(df_e18.Ideology;
                    ordered = true,
                    levels  = [-1, 0, 1, 99])


    # PT  =====================================================================================================================
    pairs = (96.0 => 99.0, 97.0 => 99.0, 98.0 => 99.0)


    replace!(df_e18.Q1513, pairs...)





    replace!(x -> ismissing(x) ? x : x < 5 ? 0.0 : x <= 10 ? 1.0 : 99.0,
                 df_e18[!, :Q1513])

    df_e18.PT = df_e18.Q1513

    return(df_e18)                                        
end