using Diversity
using EcoBase
using DataFrames
@static if VERSION < v"0.7.0-"
const deletecols! = delete!
end

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
    md[:diversity] = md[:diversity] ./
        qD(reshape(mapslices(sum, proportions, dims=(1,)),
                   size(proportions)[2]), qs)
    md[:measure] = "JostAlpha"
    return md
end

function jostalpha(asm::EcoBase.AbstractAssemblage, qs)
    hassimilarity(asm) && error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return jostalpha(occurrences(asm), qs)
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
    j = join(md, ja, on=[:q, :type_level, :type_name,
                         :partition_level, :partition_name, :div_type],
             makeunique=true)
    j[:diversity] = j[:diversity] ./ j[:diversity_1]
    j[:measure] = "JostBeta"
    deletecols!(j, [:diversity_1, :measure_1])
    return j
end

function jostbeta(asm::EcoBase.AbstractAssemblage, qs)
    hassimilarity(asm) && error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return jostbeta(occurrences(asm), qs)
end
