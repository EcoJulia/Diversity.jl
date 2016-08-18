using Diversity
using Diversity.DiversityMeasure
using Diversity.α, Diversity.ᾱ, Diversity.γ

"""
### Calculate a generalised version of richness

Calculates (species) richness of a series of columns representing
independent subcommunity counts, which is diversity at q = 0 for any
diversity measure (passed as the second argument). It also includes a
similarity matrix for the species

#### Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `DM`: diversity measure to use (one of α, ᾱ, β, β̄, ρ, ρ̄, γ)

- `proportions`: population proportions

- `Z`: similarity matrix

### Returns:
- diversity (at ecosystem level) or diversities (of subcommunities)
"""
function generalisedrichness{DM <: DiversityMeasure,
    Arr <: AbstractArray, Mat <: AbstractMatrix}(level::DiversityLevel,
                                                 ::Type{DM},
                                                 proportions::Arr,
                                                 Z::Mat =
                                                 eye(size(proportions, 1)))
    level(DM(Ecosystem(proportions, Z)), 0)
end

"""
### Calculate species richness of populations

Calculates (species) richness of a series of columns representing
independent subcommunity counts, which is diversity at q = 0

#### Arguments:
- `proportions`: population proportions

#### Returns:
- diversities of subcommunities
"""
function richness{S <: AbstractFloat}(proportions::Matrix{S})
    generalisedrichness(subcommunityDiversity, ᾱ, proportions)
end

"""
### Calculate a generalised version of Shannon entropy

Calculates Shannon entropy of a series of columns representing
independent subcommunity counts, which is log(diversity) at q = 1 for
any diversity measure (passed as the second argument). It also
includes a similarity matrix for the species

#### Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `DM`: diversity measure to use (one of α, ᾱ, β, β̄, ρ, ρ̄, γ)

- `proportions`: population proportions

- `Z`: similarity matrix

#### Returns:
- entropy (at ecosystem level) or entropies (of subcommunities)
"""
function generalisedshannon{DM <: DiversityMeasure,
    Arr <: AbstractArray, Mat <: AbstractMatrix}(level::DiversityLevel,
                                                 ::Type{DM},
                                                 proportions::Arr,
                                                 Z::Mat =
                                                 eye(size(proportions, 1)))
    log(level(DM(Ecosystem(proportions, Z)), 1))
end

"""
### Calculate Shannon entropy of populations

Calculates shannon entropy of a series of columns representing
independent subcommunity counts, which is log(diversity) at q = 1

#### Arguments:
- `proportions`: population proportions

#### Returns:
- entropies of subcommunities
"""
function shannon(proportions)
    generalisedshannon(subcommunityDiversity, ᾱ, proportions)
end

"""
### Calculate a generalised version of Simpson's index

Calculates Simpson's index of a series of columns representing
independent subcommunity counts, which is 1 / diversity at q = 2 for
any diversity measure (passed as the second argument). It also
includes a similarity matrix for the species

#### Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `DM`: diversity measure to use (one of α, ᾱ, β, β̄, ρ, ρ̄, γ)

- `proportions`: population proportions

- `Z`: similarity matrix

#### Returns:
- concentration (at ecosystem level) or concentrations (of subcommunities)
"""
function generalisedsimpson{DM <: DiversityMeasure,
    Arr <: AbstractArray, Mat <: AbstractMatrix}(level::DiversityLevel,
                                                 ::Type{DM},
                                                 proportions::Arr,
                                                 Z::Mat =
                                                 eye(size(proportions, 1)))
    level(DM(Ecosystem(proportions, Z)), 2) .^ -1
end

"""
### Calculate Simpson's index

Calculates Simpson's index of a series of columns representing
independent subcommunity counts, which is 1 / diversity (or
concentration) at q = 2

#### Arguments:
- `proportions`: population proportions

#### Returns:
- concentrations of subcommunities
"""
function simpson(proportions)
    generalisedsimpson(subcommunityDiversity, ᾱ, proportions)
end

"""
### Calculate a generalised version of the Jaccard index

Calculates a generalisation of the Jaccard index of a series of
columns representing subcommunity counts. This evaluates to is alpha / gamma
for a series of orders, repesented as a vector of qs (or a single
number).  It also includes a similarity matrix for the species. This
gives measure of the average distinctiveness of the subcommunities.

#### Arguments:
- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:
- Jaccard-related distinctivess measures
"""
function generalisedjaccard{Arr <: AbstractArray,
    Mat <: AbstractMatrix}(proportions::Arr, qs,
                           Z::Mat = eye(size(proportions, 1)))
    eco = Ecosystem(proportions, Z)
    length(eco) == 2 ||
    error("Can only calculate Jaccard index for 2 subcommunities")
    (supercommunityDiversity(α(eco), qs) ./
     supercommunityDiversity(γ(eco), qs)) - 1
end

"""
### Calculate the Jaccard index

Calculates Jaccard index (Jaccard similarity coefficient) of two
columns representing independent subcommunity counts, which is
alpha(proportions, 0) / gamma(proportions, 0) - 1

#### Arguments:
- `proportions`: population proportions

#### Returns:
- the Jaccard index
"""
function jaccard(proportions)
    generalisedjaccard(proportions, 0)[1]
end
