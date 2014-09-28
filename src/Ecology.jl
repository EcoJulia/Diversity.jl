using Diversity

## generalisedrichness() - Calculate a generalised version of richness
##
## Calculates (species) richness of a series of columns representing
## independent subcommunity counts, which is diversity at q = 0 for any
## diversity measure (passed as the second argument). It also includes
## a similarity matrix for the species
##
## Arguments:
## - measure - diversity measure to use (one of α, ᾱ, β, β̄, γ or γ̄)
## - proportions - population proportions
## - Z - similarity matrix
##
## Returns:
## - diversity (at ecosystem level) or diversities (of subcommunities)
function generalisedrichness{S <: FloatingPoint}(measure::Function,
                                                 proportions::Matrix{S},
                                                 Z::Matrix{S} =
                                                 eye(size(proportions, 1)))
    measure(proportions, 0, Z)
end

## richness() - Calculate species richness of populations
##
## Calculates (species) richness of a series of columns representing
## independent subcommunity counts, which is diversity at q = 0
##
## Arguments:
## - proportions - population proportions
##
## Returns:
## - diversities of subcommunities
function richness{S <: FloatingPoint}(proportions::Matrix{S})
    generalisedrichness(ᾱ, proportions)
end

## generalisedshannon() - Calculate a generalised version of Shannon entropy
##
## Calculates Shannon entropy of a series of columns representing
## independent subcommunity counts, which is log(diversity) at q = 1 for
## any diversity measure (passed as the second argument). It also
## includes a similarity matrix for the species
##
## Arguments:
## - measure - diversity measure to use (one of α, ᾱ, β, β̄, γ or γ̄)
## - proportions - population proportions
## - Z - similarity matrix
##
## Returns:
## - entropy (at ecosystem level) or entropies (of subcommunities)
function generalisedshannon{S <: FloatingPoint}(measure::Function,
                                                proportions::Matrix{S},
                                                Z::Matrix{S} =
                                                eye(size(proportions)[1]))
    log(measure(proportions, 1, Z))
end

## shannon() - Calculate shannon entropy of populations
##
## Calculates shannon entropy of a series of columns representing
## independent subcommunity counts, which is log(diversity) at q = 1
##
## Arguments:
## - proportions - population proportions
##
## Returns:
## - entropies of subcommunities
function shannon{S <: FloatingPoint}(proportions::Matrix{S})
    generalisedshannon(ᾱ, proportions)
end

## generalisedsimpson() - Calculate a generalised version of Simpson's index
##
## Calculates Simpson's index of a series of columns representing
## independent subcommunity counts, which is 1 / diversity at q = 2 for
## any diversity measure (passed as the second argument). It also
## includes a similarity matrix for the species
##
## Arguments:
## - measure - diversity measure to use (one of α, ᾱ, β, β̄, γ or γ̄)
## - proportions - population proportions
## - Z - similarity matrix
##
## Returns:
## - concentration (at ecosystem level) or concentrations (of subcommunities)
function generalisedsimpson{S <: FloatingPoint}(measure::Function ,
                                                proportions::Matrix{S},
                                                Z::Matrix{S} =
                                                eye(size(proportions)[1]))
    measure(proportions, 2, Z) .^ -1
end

## simpson() - Calculate Simpson's index
##
## Calculates Simpson's index of a series of columns representing
## independent subcommunity counts, which is 1 / diversity (or
## concentration) at q = 2
##
## Arguments:
## - proportions - population proportions
##
## Returns:
## - concentrations of subcommunities
function simpson{S <: FloatingPoint}(proportions::Matrix{S})
    generalisedsimpson(ᾱ, proportions)
end

## generalisedjaccard() - Calculate a generalised version of Jaccard's index
##
## Calculates a generalisation of Jaccard's index of a series of
## columns representing subcommunity counts. This evaluates to is G /
## A for a series of orders, repesented as a vector of qs (or a single
## number). It also includes a similarity matrix for the species. This
## gives measure of the average distinctiveness of the subcommunities.
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - Jaccard-related distinctivess measures
function generalisedjaccard(proportions::Matrix, qs,
                            Z::Matrix = eye(size(proportions, 1)))
    size(proportions, 2) == 2 ||
    error("Can only calculate Jaccard index for 2 subcommunities")
    A(proportions, qs, Z) ./ G(proportions, qs, Z) - 1
end

    ## jaccard() - Calculate Jaccard index
##
## Calculates Jaccard index (Jaccard similarity coefficient) of two
## columns representing independent subcommunity counts, which is
## A(proportions, 0) / G(proportions, 0) - 1
##
## Arguments:
## - proportions - population proportions
##
## Returns:
## - Jaccard index
function jaccard{S <: FloatingPoint}(proportions::Matrix{S})
    generalisedjaccard(proportions, 0)[1]
end

