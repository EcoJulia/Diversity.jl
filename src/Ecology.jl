using Diversity
using Diversity.ShortNames
using Diversity.API
using LinearAlgebra
using DataFrames
using EcoBase
using EcoBase: AbstractAssemblage

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

function richness(asm::AbstractAssemblage)
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

function shannon(asm::AbstractAssemblage)
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

function simpson(asm::AbstractAssemblage)
    hassimilarity(asm) &&
    error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return simpson(occurrences(asm))
end

"""
    generalisedjaccard(proportions::AbstractArray, qs, Z::AbstractMatrix)
    generalisedjaccard(proportions::AbstractArray, qs, sim::AbstractTypes)
    generalisedjaccard(meta::AbstractAssemblage, qs)

Calculates a generalisation of the Jaccard similarity of two columns
representing the counts of two subcommunities. This evaluates to raw
alpha / gamma - 1 for a series of orders, repesented as a vector of qs
(or a single number). It also includes an optional similarity matrix
for the species. This gives a measure of the distinctness of the
subcommunities, though we believe that beta and normalised beta have
better properties.

#### Arguments:

- `proportions`: population proportions
- `meta`: metacommunity / assemblage

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

#### Returns:

- Jaccard-related distinctivess measures
"""
function generalisedjaccard end

generalisedjaccard(proportions::AbstractMatrix, Z::AbstractMatrix) =
    generalisedjaccard(proportions, GeneralTypes(Z))

generalisedjaccard(proportions::AbstractMatrix, sim::AbstractTypes) =
    generalisedjaccard(Metacommunity(proportions, sim))

function generalisedjaccard(meta::AbstractAssemblage)
    countsubcommunities(meta) == 2 ||
    error("Can only calculate Jaccard index for 2 subcommunities")
    num = sum(minimum(getordinariness!(meta), dims = 2))
    denom = sum(maximum(getordinariness!(meta), dims = 2))
    jac = metadiv(Gamma(meta), 0)
    jac[!,:diversity] .= num / denom
    jac[!,:measure] .= "Jaccard"
    select!(jac, Not([:q]))
    return jac
end

"""
    jaccard(proportions::AbstractMatrix)
    jaccard(asm::AbstractAssemblage)

Calculates Jaccard similarity coefficient of two
columns representing independent subcommunity counts

#### Arguments:

- `proportions`: population proportions
- `asm`: assemblage / metacommunity

#### Returns:

- the Jaccard index
"""
jaccard(proportions::AbstractMatrix) = jaccard(Metacommunity(proportions))

function jaccard(asm::AbstractAssemblage)
    hassimilarity(asm) &&
    error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return generalisedjaccard(asm)
end

#=
"""
    generalisedpielou(level::DiversityLevel,
                      proportions::AbstractArray,
                      Z::AbstractMatrix)
    generalisedpielou(level::DiversityLevel,
                      proportions::AbstractArray[,
                      sim::AbstractTypes])
    generalisedpielou(level::DiversityLevel,
                      asm::AbstractAssemblage)

Calculates a generalisation of Pielou's evenness for columns
representing the counts or proportions of subcommunities. Values range from 
zero to one, with one representing complete evenness within the 
community. Since this is calculated as H / Hmax, this uses Shannon entropy
and q is effectively 1.

#### Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

#### Returns:
- Pielou's evenness metric (at metacommunity level) or metrics (of subcommunities)
"""
=#
function generalisedpielou end
generalisedpielou(level::DiversityLevel,
                  proportions::AbstractArray,
                  Z::AbstractMatrix) =
    generalisedpielou(level, Metacommunity(proportions, Z))

generalisedpielou(level::DiversityLevel,
                  proportions::AbstractArray,
                  sim::AbstractTypes = UniqueTypes(size(proportions, 1))) =
    generalisedpielou(level, Metacommunity(proportions, sim))

function generalisedpielou(level::DiversityLevel,
                           mc::AbstractAssemblage)
    hassimilarity(mc) &&
    error("Can't calculate Pielou evenness for $(typeof(gettypes(mc))) type as ill-defined maximum entropy")

    if (level == subcommunityDiversity)
        dm = ᾱ
        hmax = log.(level(dm(mc), 0).diversity)
    elseif (level == metacommunityDiversity)
        dm = Gamma
        mcab = Float64[x > 0 for x in getmetaabundance(mc)]
        mcab ./= sum(mcab)
        mc1 = Metacommunity(mcab, gettypes(mc))
        hmax = first(log.(level(dm(mc1), 0).diversity))
    else
        error("Can't calculate Pielou for $level")
    end
    gp = level(dm(mc), 1)
    gp[!,:diversity] .= log.(gp[!,:diversity]) ./ hmax
    gp[!,:measure] .= "Pielou"
    select!(gp, Not(:q))
    return gp
end

"""
    pielou(proportions::AbstractMatrix)
    pielou(asm::AbstractAssemblage)

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
    generalisedpielou(subcommunityDiversity, Metacommunity(proportions))

function pielou(asm::AbstractAssemblage)
    hassimilarity(asm) &&
    error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return generalisedpielou(subcommunityDiversity, asm)
end

"""
    gower(proportions::AbstractMatrix; countzeros::Bool = false, logscale::Bool = true)
    gower(asm::AbstractAssemblage; countzeros::Bool = false, logscale::Bool = true)

Calculates Gower's dissimarity of up to two columns representing independent subcommunity counts.
    
#### Arguments:
    
- `proportions`: population proportions; or
- `count`: population counts; or
- `asm`: Abstract Assemblage
- ``

#### Returns:
    
- Gower dissimilarity of the subcommunities
"""
function gower end

gower(proportions::AbstractArray; countzeros::Bool = false, logscale::Bool = false, normalise::Bool = countzeros) =
    gower(Metacommunity(proportions), countzeros = countzeros, logscale = logscale, normalise = normalise)

function gower(asm::AbstractAssemblage; countzeros::Bool = false, logscale::Bool = false, normalise::Bool = countzeros)
    countsubcommunities(asm) == 2 ||
    error("Can only calculate Gower distances for 2 subcommunities")

    g = meta_gamma(asm, 0)
    nz = countzeros ? nthings(asm) : noccurring(asm)
    occ = logscale ? [x == 0 ? 0 : log10(x) for x in getabundance(asm, true)] : getabundance(asm, true)
    diff = normalise ? sum(x > 0 ? 1 : 0 for x in abs.(occ[:, 1] .- occ[:, 2])) : sum(abs.(occ[:, 1] .- occ[:, 2]))
    g[!, :diversity] .= diff / nz
    g[!, :measure] .= "Gower"
    g[!, :countzeros] .= countzeros
    g[!, :logscale] .= logscale
    g[!, :normalise] .= normalise
    select!(g, Not(:q))
    return g
end
