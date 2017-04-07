using Diversity
using Diversity.ShortNames

"""
    generalisedrichness(level::DiversityLevel, proportions::AbstractArray, Z::AbstractMatrix)

Calculates species richness (diversity at q = 0) of a series of
columns representing subcommunity counts, allowing a similarity matrix
for the types / species.

# Arguments:

- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix

# Returns:

- diversity (at ecosystem level) or diversities (of subcommunities)
"""
function generalisedrichness(level::DiversityLevel,
                             proportions::AbstractArray,
                             Z::AbstractMatrix = eye(size(proportions, 1)))
    if (level == subcommunityDiversity)
        dm = ᾱ
    elseif (level == metacommunityDiversity)
        dm = Gamma
    else
        error("Can't calculate richness for $level")
    end
    gr=level(dm(Metacommunity(proportions, Z)), 0)
    gr[:measure] = "Richness"
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
    generalisedrichness(subcommunityDiversity, proportions)
end

"""
    generalisedshannon(level::DiversityLevel, proportions::AbstractArray, Z::AbstractMatrix)

Calculates Shannon entropy (log of diversity at q = 1) of a series of
columns representing independent subcommunity counts, allowing a
similarity matrix for the types / species.

# Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix

#### Returns:
- entropy (at metacommunity level) or entropies (of subcommunities)
"""
function generalisedshannon(level::DiversityLevel,
                            proportions::AbstractArray,
                            Z::AbstractMatrix = eye(size(proportions, 1)))
    if (level == subcommunityDiversity)
        dm = ᾱ
    elseif (level == metacommunityDiversity)
        dm = Gamma
    else
        error("Can't calculate richness for $level")
    end
    gs = level(dm(Metacommunity(proportions, Z)), 1)
    gs[:diversity] = log.(gs[:diversity])
    gs[:measure] = "Shannon"
    return gs
end

"""
### Calculate Shannon entropy of populations

Calculates shannon entropy (log of diversity at q = 1) of a series of
columns representing independent subcommunity counts.

#### Arguments:
- `proportions`: population proportions

#### Returns:
- entropies of subcommunities
"""
function shannon(proportions::AbstractVecOrMat)
    generalisedshannon(subcommunityDiversity, proportions)
end

"""
### Calculate a generalised version of Simpson's index

Calculates Simpson's index (1 / diversity at q = 2) of a series of
columns representing independent subcommunity counts, allowing a
similarity matrix for the types / species.

#### Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix

#### Returns:
- concentration (at ecosystem level) or concentrations (of subcommunities)
"""
function generalisedsimpson(level::DiversityLevel,
                            proportions::AbstractArray,
                            Z::AbstractMatrix = eye(size(proportions, 1)))
    if (level == subcommunityDiversity)
        dm = ᾱ
    elseif (level == metacommunityDiversity)
        dm = Gamma
    else
        error("Can't calculate richness for $level")
    end
    gs = level(dm(Metacommunity(proportions, Z)), 2)
    gs[:diversity] = gs[:diversity] .^ -1
    gs[:measure] = "Simpson"
    return gs
end

"""
    simpson(proportions::AbstractMatrix)

Calculates Simpson's index (1 / diversity at q = 2) of a series of
columns representing independent subcommunity counts.

# Arguments:

- `proportions`: population proportions

# Returns:

- concentrations of subcommunities
"""
function simpson(proportions::AbstractVecOrMat)
    generalisedsimpson(subcommunityDiversity, proportions)
end

"""
    generalisedjaccard(proportions::AbstractArray, qs, Z::AbstractMatrix)

Calculates a generalisation of the Jaccard index of two columns
representing the counts of two subcommunities. This evaluates to raw
alpha / gamma - 1 for a series of orders, repesented as a vector of qs
(or a single number). It also includes an optional similarity matrix
for the species. This gives a measure of the distinctness of the
subcommunities, though we believe that beta and normalised beta have
better properties.

# Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

# Returns:

- Jaccard-related distinctivess measures
"""
function generalisedjaccard(proportions::AbstractArray, qs,
                            Z::AbstractMatrix = eye(size(proportions, 1)))
    meta = Metacommunity(proportions, Z)
    length(meta) == 2 ||
    error("Can only calculate Jaccard index for 2 subcommunities")
    ab = metadiv(α(meta), qs)
    g = metadiv(Γ(meta), qs)
    j = join(ab, g, on=[:q, :type_level, :type_name, :partition_level, :partition_name])
    j[:diversity] = j[:diversity] ./ j[:diversity_1] - 1
    j[:measure] = "Jaccard"
    delete!(j, [:diversity_1, :measure_1])
    return j
end

"""
    jaccard(proportions::AbstractMatrix)

Calculates Jaccard index (Jaccard similarity coefficient) of two
columns representing independent subcommunity counts, which is
normmetaalpha(proportions, 0) / metagamma(proportions, 0) - 1

# Arguments:

- `proportions`: population proportions

# Returns:

- the Jaccard index
"""
function jaccard(proportions::AbstractMatrix)
    generalisedjaccard(proportions, 0)
end
