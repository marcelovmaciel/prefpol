
function ranking_to_perm(ranking)
    # Sort keys by their rank value, return ordered tuple of keys
    sorted = sort(collect(ranking); by = x -> x[2])
    return Tuple(x[1] for x in sorted)
end


function perms_out_of_rankings(profile)
    return map(ranking_to_perm, profile)
end


# TODO: check this a millionth times 

function get_consensus_ranking(profile)
    candidates = collect(keys(profile[1]))
    m = length(candidates)
    # 2. Assign each candidate an integer label
    candidate_to_label = Dict(c => i for (i, c) in enumerate(candidates))
    label_to_candidate = Dict(v => k for (k,v) in candidate_to_label)
    myperms = perms_out_of_rankings(profile)
    relabeled_perm= map(ranking->map(x-> candidate_to_label[x], ranking), myperms)
    input_for_permallows = reduce(vcat,(map(r->Vector(collect(r))', relabeled_perm)))

    @rput input_for_permallows
    R"""
    library(PerMallows)
    result <- lmm(input_for_permallows, dist.name = "cayley", estimation = "exact")
    theta <- result$theta
    mode <- result$mode
    """
    @rget mode theta

    consensus_back_from_permallows = [label_to_candidate[i] for i in mode ]
    consensus_dict = Dict(c => r for (r,c) in enumerate(consensus_back_from_permallows))
    return consensus_back_from_permallows, consensus_dict
end 





