using Diversity

"""
    jostalpha(proportions::AbstractMatrix, qs)

Calculates Jost's alpha diversity of a series of columns representing
independent community counts, for a series of orders, repesented as a
vector of qs. This is just the naive-community ecosystem diversity
divided by the naive-community beta diversity.

# Arguments:

- `proportions` relative proportions of different individuals / species
                in population (vector, or matrix where columns are
                for individual sub-communities)

- `qs` single number or vector of orders of diversity measurement

# Returns:

- DataFrame of diversities
"""
function jostalpha(proportions::AbstractMatrix, qs)
    md = metacommunityDiversity(RawAlpha(Metacommunity(proportions)), qs)
    md[:diversity] = md[:diversity] ./ qD(reshape(mapslices(sum, proportions, (1,)),
                                                  size(proportions)[2]), qs)
    md[:measure] = "JostAlpha"
    return md
end

"""
    jostbeta(proportions::AbstractMatrix, qs)

Calculates Jost's beta diversity of a series of columns representing
independent community counts, for a series of orders, repesented as a
vector of qs. This is just the naive gamma diversity divided by
Jost's alpha diversity

# Arguments:

- `proportions` relative proportions of different individuals / species
                in population (vector, or matrix where columns are
                for individual sub-communities)

- `qs` single number or vector of orders of diversity measurement

# Returns:

- DataFrame of diversities
"""
function jostbeta(proportions::AbstractMatrix, qs)
    md = metacommunityDiversity(Gamma(Metacommunity(proportions)), qs)
    ja = jostalpha(proportions, qs)
    j = join(md, ja, on=[:q, :type_level, :type_name, :partition_level, :partition_name])
    j[:diversity] = j[:diversity] ./ j[:diversity_1]
    j[:measure] = "JostBeta"
    delete!(j, [:diversity_1, :measure_1])
    return j
end
