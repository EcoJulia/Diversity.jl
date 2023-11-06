using Diversity
using Diversity.ShortNames
using Diversity.API
using LinearAlgebra
using DataFrames
using EcoBase

"""
    generalisedrichness(level::DiversityLevel, proportions::AbstractArray,
                        Z::AbstractMatrix)
    generalisedrichness(level::DiversityLevel, proportions::AbstractArray,
                        sim::AbstractTypes)

Calculates species richness (diversity at q = 0) of a series of
columns representing subcommunity counts, allowing a similarity matrix
for the types / species.

# Arguments:

- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

# Returns:

- diversity (at ecosystem level) or diversities (of subcommunities)
"""
generalisedrichness(level::DiversityLevel,
                    proportions::AbstractArray,
                    Z::AbstractMatrix = Matrix(1.0I, size(proportions, 1), size(proportions, 1))) =
    generalisedrichness(level, proportions, GeneralTypes(Z))

function generalisedrichness(level::DiversityLevel,
                             proportions::AbstractArray,
                             sim::AbstractTypes)
    if (level == subcommunityDiversity)
        dm = ᾱ
    elseif (level == metacommunityDiversity)
        dm = Gamma
    else
        error("Can't calculate richness for $level")
    end
    gr=level(dm(Metacommunity(proportions, sim)), 0)
    gr[!,:measure] .= "Richness"
    return gr
end

"""
    richness(proportions::AbstractMatrix)

Calculates species richness (diversity at q = 0) of a series of
columns representing independent subcommunity counts.

#### Arguments:
- `proportions`: population proportions

#### Returns:
- diversities of subcommunities
"""
function richness(proportions::AbstractVecOrMat)
    gr = generalisedrichness(subcommunityDiversity, proportions,
                             UniqueTypes(size(proportions, 1)))
    gr[!,:diversity] .= Int.(round.(gr[!,:diversity]))
    return gr
end

function richness(asm::EcoBase.AbstractAssemblage)
    hassimilarity(asm) && error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return richness(occurrences(asm))
end

"""
    generalisedshannon(level::DiversityLevel, proportions::AbstractArray,
                       Z::AbstractMatrix)
    generalisedshannon(level::DiversityLevel, proportions::AbstractArray,
                       sim::AbstractTypes)

Calculates Shannon entropy (log of diversity at q = 1) of a series of
columns representing independent subcommunity counts, allowing a
similarity matrix for the types / species.

# Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

#### Returns:
- entropy (at metacommunity level) or entropies (of subcommunities)
"""
function generalisedshannon end

generalisedshannon(level::DiversityLevel,
                   proportions::AbstractArray,
                   Z::AbstractMatrix = Matrix(1.0I, size(proportions, 1), size(proportions, 1))) =
    generalisedshannon(level, proportions, GeneralTypes(Z))

function generalisedshannon(level::DiversityLevel,
                            proportions::AbstractArray,
                            sim::AbstractTypes)
    if (level == subcommunityDiversity)
        dm = ᾱ
    elseif (level == metacommunityDiversity)
        dm = Gamma
    else
        error("Can't calculate richness for $level")
    end
    gs = level(dm(Metacommunity(proportions, sim)), 1)
    gs[!,:diversity] .= log.(gs[!,:diversity])
    gs[!,:measure] .= "Shannon"
    select!(gs, Not(:q))
    return gs
end

"""
    shannon(proportions::AbstractVecOrMat)

Calculates shannon entropy (log of diversity at q = 1) of a series of
columns representing independent subcommunity counts.

#### Arguments:
- `proportions`: population proportions

#### Returns:
- entropies of subcommunities
"""
shannon(proportions::AbstractVecOrMat) =
    generalisedshannon(subcommunityDiversity, proportions,
                       UniqueTypes(size(proportions, 1)))

function shannon(asm::EcoBase.AbstractAssemblage)
    hassimilarity(asm) && error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return shannon(occurrences(asm))
end

"""
    generalisedsimpson(level::DiversityLevel, proportions::AbstractArray,
                       Z::AbstractMatrix)
    generalisedsimpson(level::DiversityLevel, proportions::AbstractArray,
                       sim::AbstractTypes)

Calculates Simpson's index (1 / diversity at q = 2) of a series of
columns representing independent subcommunity counts, allowing a
similarity matrix for the types / species.

#### Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

#### Returns:
- concentration (at ecosystem level) or concentrations (of subcommunities)
"""
function generalisedsimpson end

generalisedsimpson(level::DiversityLevel,
                   proportions::AbstractArray,
                   Z::AbstractMatrix = Matrix(1.0I, size(proportions, 1), size(proportions, 1))) =
    generalisedsimpson(level, proportions, GeneralTypes(Z))

function generalisedsimpson(level::DiversityLevel,
                            proportions::AbstractArray,
                            sim::AbstractTypes)
    if (level == subcommunityDiversity)
        dm = ᾱ
    elseif (level == metacommunityDiversity)
        dm = Gamma
    else
        error("Can't calculate richness for $level")
    end
    gs = level(dm(Metacommunity(proportions, sim)), 2)
    gs[!,:diversity] .= gs[!,:diversity] .^ -1
    gs[!,:measure] .= "Simpson"
    select!(gs, Not(:q))
    return gs
end

"""
    simpson(proportions::AbstractMatrix)

Calculates Simpson's index (1 / diversity at q = 2) of a series of
columns representing independent subcommunity counts.

#### Arguments:

- `proportions`: population proportions

#### Returns:

- concentrations of subcommunities
"""
simpson(proportions::AbstractVecOrMat) =
    generalisedsimpson(subcommunityDiversity, proportions,
                       UniqueTypes(size(proportions, 1)))

function simpson(asm::EcoBase.AbstractAssemblage)
    hassimilarity(asm) &&
    error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return simpson(occurrences(asm))
end

"""
    generalisedjaccard(proportions::AbstractArray, qs, Z::AbstractMatrix)
    generalisedjaccard(proportions::AbstractArray, qs, sim::AbstractTypes)

Calculates a generalisation of the Jaccard index of two columns
representing the counts of two subcommunities. This evaluates to raw
alpha / gamma - 1 for a series of orders, repesented as a vector of qs
(or a single number). It also includes an optional similarity matrix
for the species. This gives a measure of the distinctness of the
subcommunities, though we believe that beta and normalised beta have
better properties.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

#### Returns:

- Jaccard-related distinctivess measures
"""
function generalisedjaccard end

generalisedjaccard(proportions::AbstractMatrix, qs,
                   Z::AbstractMatrix = Matrix(1.0I, size(proportions, 1),
                                              size(proportions, 1))) =
    generalisedjaccard(proportions, qs, GeneralTypes(Z))

function generalisedjaccard(proportions::AbstractMatrix, qs,
                            sim::AbstractTypes)
    meta = Metacommunity(proportions, sim)
    countsubcommunities(meta) == 2 ||
    error("Can only calculate Jaccard index for 2 subcommunities")
    ab = metadiv(α(meta), qs)
    g = metadiv(Γ(meta), qs)
    j = innerjoin(ab, g, on=[:q, :type_level, :type_name,
                             :partition_level, :partition_name, :div_type],
                  makeunique=true)
    j[!,:diversity] .= j[!,:diversity] ./ j[!,:diversity_1] .- 1
    j[!,:measure] .= "Jaccard"
    select!(j, Not([:diversity_1, :measure_1]))
    return j
end

"""
    jaccard(proportions::AbstractMatrix)

Calculates Jaccard index (Jaccard similarity coefficient) of two
columns representing independent subcommunity counts, which is
normmetaalpha(proportions, 0) / metagamma(proportions, 0) - 1

#### Arguments:

- `proportions`: population proportions

#### Returns:

- the Jaccard index
"""
jaccard(proportions::AbstractMatrix) =
    generalisedjaccard(proportions, 0,
                       UniqueTypes(size(proportions, 1)))

function jaccard(asm::EcoBase.AbstractAssemblage)
    hassimilarity(asm) &&
    error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return jaccard(occurrences(asm))
end

"""
generalisedpielou::DiversityLevel, proportions::AbstractArray,
                       Z::AbstractMatrix)
generalisedpielou(level::DiversityLevel, proportions::AbstractArray,
                       sim::AbstractTypes)

Calculates a generalisation of Pielou's evenness for columns
representing the counts or proportions of subcommunities. Values range from 
zero to one, with one representing complete evenness within the 
community. Since this is calculated as H / Hmax, and is just a proportion,
values remain unchanged regardless of the value(s) of q supplied. 

#### Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

#### Returns:
- Pielou's evenness metric (at metacommunity level) or metrics (of subcommunities)
"""
function generalisedpielou end
generalisedpielou(level::DiversityLevel,
                  proportions::AbstractArray,
                  Z::AbstractMatrix = Matrix(1.0I, size(proportions, 1),
                                             size(proportions, 1))) =
    generalisedpielou(level, proportions, GeneralTypes(Z))

function generalisedpielou(level::DiversityLevel,
                           proportions::AbstractArray,
                           sim::AbstractTypes)
    if (level == subcommunityDiversity)
        dm = ᾱ
        ns = vec(sum(x->x>0, proportions, dims=1))
    elseif (level == metacommunityDiversity)
        dm = Gamma
        ns = sum(x->x>0, proportions)
    else
        error("Can't calculate Pielou for $level")
    end
    gp = level(dm(Metacommunity(proportions, sim)), 1)
    gp[!,:diversity] .= log.(gp[!,:diversity])./log.(ns)
    gp[!,:measure] .= "Pielou"
    select!(gp, Not(:q))
    return gp
end

"""
    pielou(proportions::AbstractMatrix)

Calculates Pielou's evenness of a series of
columns representing independent subcommunity counts.

#### Arguments:

- `proportions`: population proportions

#### Returns:

- evenness of subcommunities

#### Example:
```
communitymat = [10 20 30 20 0;
                10 0 50 80 10;
                60 10 90 0 0; 
                10 10 10 10 10;
                70 70 70 70 70;
                10 0 0 90 0];

pielou(communitymat)
```
"""
pielou(proportions::AbstractVecOrMat) =
    generalisedpielou(subcommunityDiversity, proportions,
                      UniqueTypes(size(proportions, 1)))

function pielou(asm::EcoBase.AbstractAssemblage)
    hassimilarity(asm) &&
    error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return pielou(occurrences(asm))
end
