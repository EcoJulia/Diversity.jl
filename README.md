# Diversity [![Build Status](https://travis-ci.org/richardreeve/Diversity.jl.svg?branch=master)](https://travis-ci.org/richardreeve/Diversity.jl) [![Coverage Status](https://img.shields.io/coveralls/richardreeve/Diversity.jl.svg)](https://coveralls.io/r/richardreeve/Diversity.jl?branch=master)

*Diversity* is a [Julia](http://www.julialang.org) package that provides
 functionality for measuring alpha, beta and gamma diversity of
 subcommunities and ecosystems. It uses the diversity measures described
 in the arXiv paper [arXiv:1404.6520
 (q-bio.QM)](http://arxiv.org/abs/1404.6520), *How to partition
 diversity*. It also provides a series of other related and older diversity
 measures through sub-modules for compatibility. This package is still
 in alpha, and so we do not guarantee its correctness, although we are
 aware of no issues with it. Please raise an issue if you find any problems.

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
## ᾱ () - Normalised similarity-sensitive subcommunity alpha diversity.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing subcommunities, and
##   last representing values of q
function ᾱ (proportions::Matrix, qs, Z::Matrix)
subcommunityalphabar = ᾱ

## Ā () - Normalised similarity-sensitive ecosystem alpha diversity.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
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
at both the subcommunity and ecosystem level. The functions exist both with
unicode names (e.g. ᾱ ()), and with matching ascii names (e.g.
subcommunityalphabar ()):

```julia
## ᾱ () - Normalised similarity-sensitive subcommunity alpha diversity.
function ᾱ (proportions::Matrix, qs, Z::Matrix)
subcommunityalphabar = ᾱ

## Ā () - Normalised similarity-sensitive ecosystem alpha diversity.
function Ā (proportions::Matrix, qs, Z::Matrix)
ecosystemAbar = Ā

## α () - Raw similarity-sensitive subcommunity alpha diversity.
function α (proportions::Matrix, qs, Z::Matrix)
subcommunityalpha = α

## A () - Raw similarity-sensitive ecosystem alpha diversity.
function A (proportions::Matrix, qs, Z::Matrix)
ecosystemA = A

## ϵ = ρ̄  () - Normalised similarity-sensitive subcommunity beta diversity.
β̄  is retained for compatibility (= 1 / ϵ), but we believe ϵ (or ρ̄) to
be the more fundamental measure.  This is the evenness of the
subcommunity.
function ϵ (proportions::Matrix, qs, Z::Matrix)
subcommunityevenness = ϵ
subcommunityrhobar = ρ̄  = ϵ
subcommunitybetabar = β̄

## E = R̄  () - Normalised similarity-sensitive ecosystem beta diversity.
B̄  is retained for compatibility (= 1 / E), but we believe E (or R̄) to
be the more fundamental measure.  This is the average evenness of the
subcommunities.
function E (proportions::Matrix, qs, Z::Matrix)
ecosystemevenness = E
ecosystemRbar = R̄  = E
ecosystemBbar = B̄

## ρ () - Raw similarity-sensitive subcommunity beta diversity.
β is retained for compatibility (= 1 / ρ), but we believe ρ to
be the more fundamental measure.  This is the redundancy of the
subcommunity.
function ρ (proportions::Matrix, qs, Z::Matrix)
subcommunityredundancy = subcommunityrho = ρ
subcommunitybeta = β

## R () - Raw similarity-sensitive ecosystem beta diversity.
B is retained for compatibility (= 1 / R), but we believe R to
be the more fundamental measure.  This is the average redundancy of
the subcommunities.
function R (proportions::Matrix, qs, Z::Matrix)
ecosystemredundancy = ecosystemR = R
ecosystemB = B

## γ̄ () - Normalised similarity-sensitive subcommunity gamma diversity.
function γ̄ (proportions::Matrix, qs, Z::Matrix)
subcommunitygammabar = γ̄

## Ḡ () - Normalised similarity-sensitive ecosystem gamma diversity.
function Ḡ (proportions::Matrix, qs, Z::Matrix)
ecosystemGbar = Ḡ

## γ () - Raw similarity-sensitive subcommunity gamma diversity.
function γ (proportions::Matrix, qs, Z::Matrix)
subcommunitygamma = γ

## G () - Raw similarity-sensitive ecosystem gamma diversity.
function G (proportions::Matrix, qs, Z::Matrix)
ecosystemG = G
```

We also provide a general function for extract any diversity measure for a
series of subcommunity relative abundances:

```julia
## diversity () - calculates subcommunity and ecosystem diversities
##
## Calculates any diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
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
##                     subcommunity diversities
## - returnweights   - boolean describing whether to return subcommunity weights
##
## Returns:
## - some or all (as tuple) of:
##   - vector of ecosystem diversities representing values of q
##   - array of diversities, first dimension representing subcommunities, and
##     last representing values of q
##   - vector of subcommunity weights
function diversity (measure::Function, proportions::Matrix, qs, Z::Matrix,
                    returnecosystem::Bool,
                    returncommunity::Bool,
                    returnweights::Bool)
```

And we can calculate the proportions that subcommunities each
contribute to ecosystem diversity per subcommunity or per individual:

```julia
## contributions () - Calculate diversity contributions from subcommunities
##
## Calculates proportions that subcommunities each contribute to
## ecosystem diversity per subcommunity (perindividual = false), or
## per individual (perindividual = true) - in the latter case scaled
## so that the total # of individuals is 1, since we only have
## relative abundances.
##
## Arguments:
## - measure - diversity measure to use - one of α , ᾱ , β , β̄ , γ or γ̄
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - perindividual - do we measure per individual in population (true)
##                   or per subcommunity (false)
## - Z - similarity matrix
## - returnecosystem - boolean describing whether to return the
##                     ecosystem diversity
## - returncommunity - boolean describing whether to return the
##                     subcommunity diversities
## - returnweights   - boolean describing whether to return subcommunity weights
##
## Returns:
## - contributions of subcommunities to ecosystem diversity (of type measure)
## - and none, some or all (in a tuple) of:
##   - vector of ecosystem diversities representing values of q
##   - array of diversities, first dimension representing subcommunities, and
##     last representing values of q
##   - vector of subcommunity weights
function contributions (measure::Function, proportions::Matrix, qs,
                        perindividual::Bool, Z::Matrix,
                        returnecosystem::Bool,
                        returncommunity::Bool,
                        returnweights::Bool)
```

The package also provides sub-modules with other diversity measures.
Old ecological diversity measures are found in the .Ecology sub-module:

```julia
using Diversity.Ecology

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
function richness(proportions)

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
function shannon(proportions)

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
function simpson(proportions)

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
function jaccard(proportions::Matrix)
```

We have also developed generalised version of these that relate to our
general measures of alpha, beta and gamma diversity at subcommunity
and ecosystem measures. The generalisations of the richness, Shannon
and Simpson are the only standard measures we are aware of whose
subcommunity components sum directly to the corresponding ecosystem
measure (although note that Simpson's index decreases for increased
diversity, so small components are more diverse):

```julia
## generalisedrichness () - Calculate a generalised version of richness
##
## Calculates (species) richness of a series of columns representing
## independent subcommunity counts, which is diversity at q = 0 for any
## diversity measure (passed as the second argument). It also includes
## a similarity matrix for the species
##
## Arguments:
## - measure - diversity measure to use - one of α , ᾱ , β , β̄ , γ or γ̄ 
## - proportions - population proportions
## - Z - similarity matrix
##
## Returns:
## - diversity (at ecosystem level) or diversities (of subcommunities)
function generalisedrichness (measure::Function,
                              proportions::Matrix,
                              Z::Matrix)

## generalisedshannon () - Calculate a generalised version of Shannon entropy
##
## Calculates Shannon entropy of a series of columns representing
## independent subcommunity counts, which is log(diversity) at q = 1 for
## any diversity measure (passed as the second argument). It also
## includes a similarity matrix for the species
##
## Arguments:
## - measure - diversity measure to use - one of α , ᾱ , β , β̄ , γ or γ̄ 
## - proportions - population proportions
## - Z - similarity matrix
##
## Returns:
## - entropy (at ecosystem level) or entropies (of subcommunities)
function generalisedshannon (measure::Function,
                             proportions::Matrix,
                             Z::Matrix)

## generalisedsimpson () - Calculate a generalised version of Simpson's index
##
## Calculates Simpson's index of a series of columns representing
## independent subcommunity counts, which is 1 / diversity at q = 2 for
## any diversity measure (passed as the second argument). It also
## includes a similarity matrix for the species
##
## Arguments:
## - measure - diversity measure to use - one of α , ᾱ , β , β̄ , γ or γ̄ 
## - proportions - population proportions
## - Z - similarity matrix
##
## Returns:
## - concentration (at ecosystem level) or concentrations (of subcommunities)
function generalisedsimpson (measure::Function,
                             proportions::Matrix,
                             Z::Matrix)

## generalisedjaccard () - Calculate a generalised version of Jaccard's index
##
## Calculates a generalisation of Jaccard's index of two columns
## representing subcommunity counts. This evaluates to is A / G - 1
## for a series of orders, repesented as a vector of qs (or a single
## number). It also includes a similarity matrix for the species. This
## gives measure of the average distinctiveness of the
## subcommunities.
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - Jaccard-related distinctivess measures
function generalisedjaccard (proportions::Matrix, qs,
                             Z::Matrix)
```

[Hill numbers](http://www.jstor.org/stable/1934352) are found in the
.Hill sub-module:

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

And Jost's
[diversity](http://dx.doi.org/10.1111/j.2006.0030-1299.14714.x)
[measures](http://www.esajournals.org/doi/abs/10.1890/06-1736.1) are
found in the .Jost sub-module:

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

## jostbeta () - calculate Jost's beta diversity of multiple subcommunities
##
## Calculates Jost's beta diversity of a series of columns
## representing independent subcommunity counts, for a series of orders,
## repesented as a vector of qs. This is just the naive gamma
## diversity divided by the naive alpha diversity
##
## Arguments:
## - proportions - relative proportions of different individuals / species
##                 in population (vector, or matrix where columns are
##                 for individual subcommunities)
## - qs - single number or vector of orders of diversity measurement
##
## Returns:
## - array of diversities, first dimension representing subcommunities, and
##   last representing values of q
function jostbeta (proportions, qs)
jostβ = jostbeta
```

