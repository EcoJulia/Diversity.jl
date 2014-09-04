# Diversity [![Build Status](https://travis-ci.org/richardreeve/Diversity.jl.svg?branch=master)](https://travis-ci.org/richardreeve/Diversity.jl) [![Coverage Status](https://img.shields.io/coveralls/richardreeve/Diversity.jl.svg)](https://coveralls.io/r/richardreeve/Diversity.jl?branch=master)

*Diversity* is a [Julia](http://www.julialang.org) package that provides
 functionality for measuring alpha, beta and gamma diversity of
 communities and ecosystems. It uses the diversity measures described
 in the arXiv paper [arXiv:1404.6520
 (q-bio.QM)](http://arxiv.org/abs/1404.6520), *How to partition
 diversity*. It also provides a series of other related and older diversity
 measures through sub-modules for compatibility.

## Install

*Diversity* is in `METADATA` and can be installed via `Pkg.add("Diversity")`.

## Usage

Accessing the main functionality in the package is simple:

```julia
using Diversity
...
diversities = ᾱ(proportions, [0, 1, 2, Inf], Z)
```

The main package provides basic diversity measures (from
[Hill, 1973](http://www.jstor.org/stable/1934352)):

```julia
## qD () - calculate Hill number / naive diversity of order q of
## population(s) with given relative proportions
##
## Arguments:
## - proportions - relative proportions of different individuals / species
##                 in population (vector, or matrix where columns are
##                 individual populations) 
## - qs - single number or vector of orders of diversity measurement
##
## Returns:
## - Diversity of order qs (single number or vector of diversities)
function qD (proportions, qs)
```

And measures which can account for similarity between individuals and
species (from [Leinster and Cobbold,
2012](http://www.esajournals.org/doi/abs/10.1890/10-2402.1)):

```julia
## qDZ () - calculates Leinster-Cobbold general diversity of order q(s)
## of population(s) with given relative proportions, and similarity
## matrix Z
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population (vector, or matrix where columns are
##                 individual populations)
## - qs - single number or vector of orders of diversity measurement
## - Z - similarity matrix
##
## Returns:
## - Diversity of order qs (single number or vector of diversities)
function qDZ (proportions, qs, Z)

```

It also provides generalised alpha, beta and gamma diversity measures at the
level of the ecosystem and its component subcommunities (from [Reeve et al,
2014](http://arxiv.org/abs/1404.6520)). The normalised alpha diversities are
described here:

```julia
## ᾱ () - Normalised similarity-sensitive sub-community alpha diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function ᾱ (proportions::Matrix, qs, Z::Matrix)
communityalphabar = ᾱ

## Ā () - Normalised similarity-sensitive ecosystem alpha diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of ecosystem diversities representing values of q
function Ā (proportions::Matrix, qs, Z::Matrix)
ecosystemAbar = Ā
```

There are also matching normalised and raw, alpha, beta and gamma diversities
at both the sub-community and ecosystem level. The functions exist both with
unicode names (e.g. ᾱ ()), and with matching ascii names (e.g.
communityalphabar ()):

```julia
## ᾱ () - Normalised similarity-sensitive sub-community alpha diversity.
function ᾱ (proportions::Matrix, qs, Z::Matrix)
communityalphabar = ᾱ

## Ā () - Normalised similarity-sensitive ecosystem alpha diversity.
function Ā (proportions::Matrix, qs, Z::Matrix)
ecosystemAbar = Ā

## α () - Raw similarity-sensitive sub-community alpha diversity.
function α (proportions::Matrix, qs, Z::Matrix)
communityalpha = α

## A () - Raw similarity-sensitive ecosystem alpha diversity.
function A (proportions::Matrix, qs, Z::Matrix)
ecosystemA = A

## β̄ () - Normalised similarity-sensitive sub-community beta diversity.
function β̄ (proportions::Matrix, qs, Z::Matrix)
communitybetabar = β̄

## B̄ () - Normalised similarity-sensitive ecosystem beta diversity.
function B̄ (proportions::Matrix, qs, Z::Matrix)
ecosystemBbar = B̄

## β () - Raw similarity-sensitive sub-community beta diversity.
function β (proportions::Matrix, qs, Z::Matrix)
communitybeta = β

## B () - Raw similarity-sensitive ecosystem beta diversity.
function B (proportions::Matrix, qs, Z::Matrix)
ecosystemB = B

## γ̄ () - Normalised similarity-sensitive sub-community gamma diversity.
function γ̄ (proportions::Matrix, qs, Z::Matrix)
communitygammabar = γ̄

## Ḡ () - Normalised similarity-sensitive ecosystem gamma diversity.
function Ḡ (proportions::Matrix, qs, Z::Matrix)
ecosystemGbar = Ḡ

## γ () - Raw similarity-sensitive sub-community gamma diversity.
function γ (proportions::Matrix, qs, Z::Matrix)
communitygamma = γ

## G () - Raw similarity-sensitive ecosystem gamma diversity.
function G (proportions::Matrix, qs, Z::Matrix)
ecosystemG = G
```

We also provide a general function for extract any diversity measure for a
series of sub-community relative abundances:

```julia
## diversity () - calculates sub-community and ecosystem diversities
##
## Calculates any diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs, with similarity matrix Z, by default the (naïve) identity
## matrix.
##
## Arguments:
## - measure - the diversity to be used - one of α , ᾱ , β , β̄ , γ or γ̄
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
## - returnecosystem - boolean describing whether to return the
##                     ecosystem diversity
## - returncommunity - boolean describing whether to return the
##                     community diversities
## - returnweights   - boolean describing whether to return community weights
##
## Returns:
## - some or all (as tuple) of:
##   - vector of ecosystem diversities representing values of q
##   - array of diversities, first dimension representing sub-communities, and
##     last representing values of q
##   - vector of community weights
function diversity (measure::Function, proportions::Matrix, qs, Z::Matrix,
                    returnecosystem::Bool,
                    returncommunity::Bool,
                    returnweights::Bool)
```

And we can calculate the proportions that sub-communities each
contribute to ecosystem diversity per sub-community or per individual:

```julia
## contributions () - Calculate diversity contributions from sub-communities
##
## Calculates proportions that sub-communities each contribute to
## ecosystem diversity per sub-community (perindividual = false), or
## per individual (perindividual = true) - in the latter case scaled
## so that the total # of individuals is 1, since we only have
## relative abundances.
##
## Arguments:
## - measure - diversity measure to use
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - perindividual - do we measure per individual in population (true)
##                   or per sub-community (false)
## - Z - similarity matrix
## - returnecosystem - boolean describing whether to return the
##                     ecosystem diversity
## - returncommunity - boolean describing whether to return the
##                     community diversities
## - returnweights   - boolean describing whether to return community weights
##
## Returns:
## - contributions of sub-communities to ecosystem diversity (of type measure)
## - and none, some or all (in a tuple) of:
##   - vector of ecosystem diversities representing values of q
##   - array of diversities, first dimension representing sub-communities, and
##     last representing values of q
##   - vector of community weights
function contributions (measure::Function, proportions::Matrix, qs,
                        perindividual::Bool, Z::Matrix,
                        returnecosystem::Bool,
                        returncommunity::Bool,
                        returnweights::Bool)
```

The package also provides sub-modules with other diversity measures.
Old diversity measures are found in the .Compatibility sub-module:

```julia
using Diversity.Compatibility

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
function richness(proportions)

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
function shannon(proportions)

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
function simpson(proportions)
```

We have also developed generalised version of these measures that
relate to our general measures of alpha, beta and gamma diversity at
sub-community and ecosystem measures. They are the only standard
measures whose sub-community components sum directly to the
corresponding ecosystem measure (although note that Simpson's index
decreases for increased diversity, so small components are more
diverse):

```julia
## generalisedrichness () - Calculate a generalised version of richness
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
function generalisedrichness (proportions::Matrix,
                              measure::Function = ᾱ,
                              Z::Matrix)

## generalisedshannon () - Calculate a generalised version of Shannon entropy
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
function generalisedshannon (proportions::Matrix,
                             measure::Function = ᾱ,
                             Z::Matrix)

## generalisedsimpson () - Calculate a generalised version of Simpson's index
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
function generalisedsimpson (proportions::Matrix,
                             measure::Function = ᾱ,
                             Z::Matrix)
```

Hill numbers are found in the .Hill sub-module:

```julia
using Diversity.Hill

## hillnumber () - calculate Hill number / naive diversity of order q of
## population(s) with given relative proportions
##
## Arguments:
## - proportions - relative proportions of different individuals / species
##                 in population (vector, or matrix where columns are
##                 individual populations) 
## - qs - single number or vector of orders of diversity measurement
##
## Returns:
## - Diversity of order qs (single number or vector of diversities)
function hillnumber (proportions, qs)
```

And Jost's diversity measures are found in the .Jost sub-module:

```julia
using Diversity.Jost

## jostD () - calculate Hill number / naive diversity of order q of
## population(s) with given relative proportions
##
## Arguments:
## - proportions - relative proportions of different individuals / species
##                 in population (vector, or matrix where columns are
##                 individual populations) 
## - qs - single number or vector of orders of diversity measurement
##
## Returns:
## - Diversity of order qs (single number or vector of diversities)
function jostD (proportions, qs)

## jostbeta () - calculate Jost's beta diversity of multiple sub-communities
##
## Calculates Jost's beta diversity of a series of columns
## representing independent community counts, for a series of orders,
## repesented as a vector of qs. This is just the naive gamma
## diversity divided by the naive alpha diversity
##
## Arguments:
## - proportions - relative proportions of different individuals / species
##                 in population (vector, or matrix where columns are
##                 for individual sub-communities)
## - qs - single number or vector of orders of diversity measurement
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function jostbeta (proportions, qs)
```

