using Diversity

## generalisedrichness() - Calculate a generalised version of richness
##
## Calculates (species) richness of a series of columns representing
## independent community counts, which is diversity at q = 0 for any
## diversity measure (passed as the second argument). It also includes
## a similarity matrix for the species
##
## Arguments:
## - proportions - population proportions
## - measure - diversity measure to use, by default ᾱ
## - Z - similarity matrix
##
## Returns:
## - diversity (at ecosystem level) or diversities (of sub-communities)
function generalisedrichness{S <: FloatingPoint}(proportions::Matrix{S},
                                                 measure::Function = ᾱ,
                                                 Z::Matrix{S} = eye(size(proportions)[1]))
    measure(proportions, 0, Z)
end

## richness() - Calculate species richness of populations
##
## Calculates (species) richness of a series of columns representing
## independent community counts, which is diversity at q = 0
##
## Arguments:
## - proportions - population proportions
##
## Returns:
## - diversities of sub-communities
function richness{S <: FloatingPoint}(proportions::Matrix{S})
    generalisedrichness(proportions, ᾱ)
end

## generalisedshannon() - Calculate a generalised version of Shannon entropy
##
## Calculates Shannon entropy of a series of columns representing
## independent community counts, which is log(diversity) at q = 1 for
## any diversity measure (passed as the second argument). It also
## includes a similarity matrix for the species
##
## Arguments:
## - proportions - population proportions
## - measure - diversity measure to use, by default ᾱ
## - Z - similarity matrix
##
## Returns:
## - entropy (at ecosystem level) or entropies (of sub-communities)
function generalisedshannon{S <: FloatingPoint}(proportions::Matrix{S},
                                                measure::Function = ᾱ,
                                                Z::Matrix{S} = eye(size(proportions)[1]))
    log(measure(proportions, 1, Z))
end

## shannon() - Calculate shannon entropy of populations
##
## Calculates shannon entropy of a series of columns representing
## independent community counts, which is log(diversity) at q = 1
##
## Arguments:
## - proportions - population proportions
##
## Returns:
## - entropies of sub-communities
function shannon{S <: FloatingPoint}(proportions::Matrix{S})
    generalisedshannon(proportions, ᾱ)
end

## generalisedsimpson() - Calculate a generalised version of Simpson's index
##
## Calculates Simpson's index of a series of columns representing
## independent community counts, which is 1 / diversity at q = 2 for
## any diversity measure (passed as the second argument). It also
## includes a similarity matrix for the species
##
## Arguments:
## - proportions - population proportions
## - measure - diversity measure to use, by default ᾱ
## - Z - similarity matrix
##
## Returns:
## - concentration (at ecosystem level) or concentrations (of sub-communities)
function generalisedsimpson{S <: FloatingPoint}(proportions::Matrix{S},
                                                measure::Function = ᾱ,
                                                Z::Matrix{S} = eye(size(proportions)[1]))
    measure(proportions, 2, Z) .^ -1
end

## simpson() - Calculate Simpson's index
##
## Calculates Simpson's index of a series of columns representing
## independent community counts, which is 1 / diversity (or
## concentration) at q = 2
##
## Arguments:
## - proportions - population proportions
##
## Returns:
## - concentrations of sub-communities
function simpson{S <: FloatingPoint}(proportions::Matrix{S})
    generalisedsimpson(proportions, ᾱ)
end

