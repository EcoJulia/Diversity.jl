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

The package also provides sub-modules with other diversity measures:

```julia
using Diversity.Compatibility
...
simp = simpson(proportions)
rich = richness(proportions)
shan = shannon(proportions)

using Diversity.Hill
...
diversities = hillnumber(proportions, [0, 1, 2, Inf])

using Diversity.Jost
...
diversities = jostbeta(proportions, [0, 1, 2, Inf])
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
## diversity() - calculates sub-community and ecosystem diversities
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
## - returnweights   - boolean describing whether to return the
##                     community weights
##
## Returns:
## - any or all (as tuple) of:
##   - vector of ecosystem diversities representing values of q
##   - array of diversities, first dimension representing sub-communities, and
##     last representing values of q
##   - vector of community weights
function diversity (measure::Function, proportions::Matrix, qs, Z::Matrix,
                    returnecosystem::Bool,
                    returncommunity::Bool,
                    returnweights::Bool)
```