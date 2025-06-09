function load_and_prepare_scores_df(data_path::String; CANDIDATOS = CANDIDATOS_eseb2022)
    # Load SPSS file
    df_e22 = load_spss_file(data_path)

    # Metadata
    PARTIDOS = ["PDT", "PL", "PODEMOS", "PP", "PT", "PSB", "PSD", "PSDB", "PSOL", "REDE", "REP", "UB", "MDB"]

    # Helper function
    build_column_symbols(base::String, n::Int) = [Symbol(base * string(i)) for i in 1:n]

    # Copy and rename
    scores_df = copy(df_e22)
    rename!(scores_df, Dict(zip(build_column_symbols("Q16_", 13), PARTIDOS)))
    rename!(scores_df, Dict(zip(build_column_symbols("Q17_", 13), CANDIDATOS)))

    # Recode D10 (religion)
    replace!(x -> x in (99.0, 100.0, 101.0, 102.0) ? 95.0 : x, scores_df.D10)
    replace!(x -> x in (97.0, 98.0) ? 97.0 : x, scores_df.D10)
    scores_df.D10 = categorical(scores_df.D10)
    scores_df.Religion = scores_df.D10

    # Recode D02 (sex)
    scores_df.D02 = categorical(scores_df.D02)
    scores_df.Sex = scores_df.D02
    # Recode D12a (race)
    replace!(x -> x in (97.0, 98.0) ? 97.0 : x, scores_df.D12a)
    
    scores_df.D12a = categorical(scores_df.D12a)
    
    scores_df.Race = scores_df.D12a

    # Recode Q19 problematic codes
    replace!(x -> x in (95.0, 96.0, 98.0) ? 99.0 : x, scores_df.Q19)

    # Recode Q19 into ideological categories
    new = Vector{Int64}(undef, nrow(scores_df))
    @inbounds for (i, x) in pairs(scores_df.Q19)
        if x == 99
            new[i] = 99
        elseif 0 ≤ x ≤ 3
            new[i] = -1
        elseif x ≤ 6
            new[i] = 0
        elseif x ≤ 10
            new[i] = 1
        else
            new[i] = 99
        end
    end
    scores_df.Ideology= categorical(new;
        ordered = true,
        levels  = [-1, 0, 1, 99])
    
    scores_df[!, :PT] = Float64[
        coalesce((x ≥ 5) && (x ≤ 10), false) ? 1.0 : 0.0
        for x in scores_df[!, :PT]
            ]
#=     pt_clean = coalesce.(scores_df.PT, 0)      # replaces missing with 0
    code = ifelse.(pt_clean .< 5,            0,
        ifelse.(pt_clean .<= 10,         1,
                                        99)) =#
   transform!(scores_df,
           :Q31_7 => ByRow(x -> x in (97.0, 98.0) ? 3.0 : x) => :Abortion)
    return scores_df
end



function load_and_prepare_e2006(df_path; candidates = candidates2006)
    df_e06 = load_spss_file(df_path)
    letters = ['a','b','c','d','e','f']

    function build_letter_column_symbols(base::AbstractString, letters::Vector{Char})
        return [Symbol(base * string(c)) for c in letters]
    end

    rename!(df_e06, Dict(zip(build_letter_column_symbols("eseb16", letters), candidates)))



    df_e06.peso = df_e06.peso_1

    df_e06.Sex = categorical(df_e06.SEXO)


    # df_e06.eseb8 |> proportionmap again, gigantic missing here 

    pt_clean = coalesce.(df_e06.eseb15a, 0)      # replaces missing with 0

    code = ifelse.(pt_clean .< 5,            0,
                ifelse.(pt_clean .<= 10,         1,
                                                99))
    df_e06.PT = code 



        # Recode Q19 into ideological categories
    new = Vector{Int64}(undef, nrow(df_e06))
            @inbounds for (i, x) in pairs(df_e06.eseb19)
                if x == 98.
                    new[i] = 98
                elseif 0 ≤ x ≤ 3
                    new[i] = -1
                elseif x ≤ 6
                    new[i] = 0
                elseif x ≤ 10
                    new[i] = 1
                else
                    new[i] = 98
                end
        end


    df_e06.Ideology= categorical(new;
                ordered = true,
                levels  = [-1, 0, 1, 98])
    return(df_e06)
end



function load_and_prepare_e2018(df_path; CANDIDATOS = candidates2018)
   

    df_e18 = load_spss_file(df_path)

    function build_column_symbols(base::AbstractString, n::Integer;
                                minwidth::Int = 2)
        width = max(minwidth, length(string(n)))          # e.g. n = 21  → width = 2
        return [Symbol(base * @sprintf("%0*d", width, i)) for i in 1:n]
    end



    rename!(df_e18, Dict(zip(build_column_symbols("Q16", 21), CANDIDATOS)))


    replace!(x -> x in (97., 98., 99.) ? 99.0 : x, df_e18.D10)
    df_e18.D10 = categorical(df_e18.D10)
    df_e18.Religion = df_e18.D10

    df_e18.D2_SEXO = categorical(df_e18.D2_SEXO)
    df_e18.Sex  = df_e18.D2_SEXO


    replace!(x -> x in (8., 9.) ? 8.0 : x, df_e18.D12A)

    df_e18.Race = df_e18.D12A



    replace!(x -> x in (95.0, 97.0, 98.0) ? 98.0 : x, df_e18.Q18)




    # Recode Q19 into ideological categories
    new = Vector{Int64}(undef, nrow(df_e18))
        @inbounds for (i, x) in pairs(df_e18.Q18)
            if x == 98.
                new[i] = 98
            elseif 0 ≤ x ≤ 3
                new[i] = -1
            elseif x ≤ 6
                new[i] = 0
            elseif x ≤ 10
                new[i] = 1
            else
                new[i] = 98
            end
    end

    df_e18.Ideology= categorical(new;
            ordered = true,
            levels  = [-1, 0, 1, 98])
        


    pt_clean = coalesce.(df_e18.Q1513, 0)      # replaces missing with 0

    code = ifelse.(pt_clean .< 5,            0,
            ifelse.(pt_clean .<= 10,         1,
                                            99))
    df_e18.PT = code
    return(df_e18)                                        
end