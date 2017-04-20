__precompile__()

"""
The main **Diversity** module provides basic numbers-equivalent
diversity measures (described in
[Hill, 1973](http://www.jstor.org/stable/1934352)),
similarity-sensitive diversity measures (generalised from Hill, and
described in
[Leinster and Cobbold, 2012](http://www.esajournals.org/doi/abs/10.1890/10-2402.1)),
and related alpha, beta and gamma diversity measures at the level of
the metacommunity and its component subcommunities (generalised in
turn from Leinster and Cobbold, and described in
[Reeve et al, 2014](http://arxiv.org/abs/1404.6520)). The diversity
functions exist both with unicode names (e.g. ```ᾱ()```), which are
not automatically exported (as we feel they are too short) and with
matching longer ASCII names (e.g. ```NormalisedAlpha()```), which are.
We also provide functions to calculate appropriate
```subcommunityDiversity()``` and ```metacommunityDiversity()```
values for each measure, a general ```diversity()``` function for
extract any diversity measure at a series of scales.
"""
module Diversity

"""
The Diversity.API submodule should be `import`ed if you want to create a
new type, partition or metacommunity subtype. Otherwise it can be
ignored.
"""
module API
include("API.jl")
export AbstractTypes, AbstractPartition, AbstractMetacommunity
export _gettypes, _getpartition
export _counttypes, _countsubcommunities, _getnames
export _getabundance, _getmetaabundance, _getweight
export _getordinariness!, _getmetaordinariness!
export _calcabundance, _calcsimilarity, _calcordinariness
export floattypes, typematch, mcmatch
export sumovertypes, sumoversubcommunities
end

include("Interface.jl")
export gettypes, getpartition
export counttypes, countsubcommunities
export gettypenames, getsubcommunitynames, getnames
export getabundance, getmetaabundance, getweight
export getordinariness!, getmetaordinariness!
export calcabundance, calcsimilarity, calcordinariness

include("Metacommunity.jl")
export Subcommunities, Onecommunity
export GeneralTypes, UniqueTypes, Species, Taxonomy
export Metacommunity
export inddiv, subdiv, metadiv

include("EffectiveNumbers.jl")
export qD, qDZ

include("DiversityMeasure.jl")
export DiversityLevel
export individualDiversity, subcommunityDiversity, metacommunityDiversity
export DiversityMeasure, PowerMeanMeasure, RelativeEntropyMeasure
export RawAlpha, NormalisedAlpha
export RawBeta, NormalisedBeta, RawRho, NormalisedRho
export Gamma

## We do not directly export ᾱ, α, β̄, β, ρ̄, ρ, γ as they're too short
## only via Diversity.ShortNames
module ShortNames
using Diversity

const α = RawAlpha
const ᾱ = NormalisedAlpha
const β = RawBeta
const β̄ = NormalisedBeta
const ρ = RawRho
const ρ̄ = NormalisedRho
const γ = Gamma

export α, ᾱ, β, β̄, ρ, ρ̄
# γ actually can't be exported like this - it'll always just be Shortnames.γ, so we export Γ instead
const Γ = Gamma
export Γ

end

export getName, getASCIIName, getFullName

include("GeneralisedDiversities.jl")
export diversity
export norm_sub_alpha, raw_sub_alpha, norm_sub_beta, raw_sub_beta
export norm_sub_rho, raw_sub_rho, sub_gamma
export norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta
export norm_meta_rho, raw_meta_rho, meta_gamma

## Deprecate short names as too ambiguous
@deprecate(Dᾱ, subcommunityalphabar)
@deprecate(Dα, subcommunityalpha)
@deprecate(Dβ̄, subcommunitybetabar)
@deprecate(Dβ, subcommunitybeta)
@deprecate(Dρ̄, subcommunityrhobar)
@deprecate(Dρ, subcommunityrho)
@deprecate(Dγ̄, subcommunitygammabar)
@deprecate(Dγ, subcommunitygamma)
@deprecate(DĀ, supercommunityAbar)
@deprecate(DA, supercommunityA)
@deprecate(DB̄, supercommunityBbar)
@deprecate(DB, supercommunityB)
@deprecate(DR̄, supercommunityRbar)
@deprecate(DR, supercommunityR)
@deprecate(DḠ, supercommunityGbar)
@deprecate(DG, supercommunityG)

## Deprecate ecosystem-related names
@deprecate(ecosystemAbar(pop::AbstractArray, qs, Z::AbstractMatrix),
           norm_meta_alpha(Metacommunity(pop, Z), qs))
@deprecate(ecosystemA(pop::AbstractArray, qs, Z::AbstractMatrix),
           raw_meta_alpha(Metacommunity(pop, Z), qs))
@deprecate(ecosystemBbar(pop::AbstractArray, qs, Z::AbstractMatrix),
           norm_meta_beta(Metacommunity(pop, Z), qs))
@deprecate(ecosystemB(pop::AbstractArray, qs, Z::AbstractMatrix),
           raw_meta_beta(Metacommunity(pop, Z), qs))
@deprecate(ecosystemRbar(pop::AbstractArray, qs, Z::AbstractMatrix),
           norm_meta_rho(Metacommunity(pop, Z), qs))
@deprecate(ecosystemR(pop::AbstractArray, qs, Z::AbstractMatrix),
           raw_meta_rho(Metacommunity(pop, Z), qs))
@deprecate(ecosystemGbar(pop::AbstractArray, qs, Z::AbstractMatrix),
           meta_gamma(Metacommunity(pop, Z), qs))
@deprecate(ecosystemG(pop::AbstractArray, qs, Z::AbstractMatrix),
           meta_gamma(Metacommunity(pop, Z), qs))

@deprecate(ecosystemAbar(pop::AbstractArray, qs),
           norm_meta_alpha(Metacommunity(pop), qs))
@deprecate(ecosystemA(pop::AbstractArray, qs),
           raw_meta_alpha(Metacommunity(pop), qs))
@deprecate(ecosystemBbar(pop::AbstractArray, qs),
           norm_meta_beta(Metacommunity(pop), qs))
@deprecate(ecosystemB(pop::AbstractArray, qs),
           raw_meta_beta(Metacommunity(pop), qs))
@deprecate(ecosystemRbar(pop::AbstractArray, qs),
           norm_meta_rho(Metacommunity(pop), qs))
@deprecate(ecosystemR(pop::AbstractArray, qs),
           raw_meta_rho(Metacommunity(pop), qs))
@deprecate(ecosystemGbar(pop::AbstractArray, qs),
           meta_gamma(Metacommunity(pop), qs))
@deprecate(ecosystemG(pop::AbstractArray, qs),
           meta_gamma(Metacommunity(pop), qs))

## Deprecate supercommunity-related names
@deprecate(supercommunityAbar(pop::AbstractArray, qs, Z::AbstractMatrix),
           norm_meta_alpha(Metacommunity(pop, Z), qs))
@deprecate(supercommunityA(pop::AbstractArray, qs, Z::AbstractMatrix),
           raw_meta_alpha(Metacommunity(pop, Z), qs))
@deprecate(supercommunityBbar(pop::AbstractArray, qs, Z::AbstractMatrix),
           norm_meta_beta(Metacommunity(pop, Z), qs))
@deprecate(supercommunityB(pop::AbstractArray, qs, Z::AbstractMatrix),
           raw_meta_beta(Metacommunity(pop, Z), qs))
@deprecate(supercommunityRbar(pop::AbstractArray, qs, Z::AbstractMatrix),
           norm_meta_rho(Metacommunity(pop, Z), qs))
@deprecate(supercommunityR(pop::AbstractArray, qs, Z::AbstractMatrix),
           raw_meta_rho(Metacommunity(pop, Z), qs))
@deprecate(supercommunityGbar(pop::AbstractArray, qs, Z::AbstractMatrix),
           meta_gamma(Metacommunity(pop, Z), qs))
@deprecate(supercommunityG(pop::AbstractArray, qs, Z::AbstractMatrix),
           meta_gamma(Metacommunity(pop, Z), qs))

@deprecate(supercommunityAbar(pop::AbstractArray, qs),
           norm_meta_alpha(Metacommunity(pop), qs))
@deprecate(supercommunityA(pop::AbstractArray, qs),
           raw_meta_alpha(Metacommunity(pop), qs))
@deprecate(supercommunityBbar(pop::AbstractArray, qs),
           norm_meta_beta(Metacommunity(pop), qs))
@deprecate(supercommunityB(pop::AbstractArray, qs),
           raw_meta_beta(Metacommunity(pop), qs))
@deprecate(supercommunityRbar(pop::AbstractArray, qs),
           norm_meta_rho(Metacommunity(pop), qs))
@deprecate(supercommunityR(pop::AbstractArray, qs),
           raw_meta_rho(Metacommunity(pop), qs))
@deprecate(supercommunityGbar(pop::AbstractArray, qs),
           meta_gamma(Metacommunity(pop), qs))
@deprecate(supercommunityG(pop::AbstractArray, qs),
           meta_gamma(Metacommunity(pop), qs))

## Deprecate subcommunity-related names
@deprecate(subcommunityalphabar(pop::AbstractArray, qs, Z::AbstractMatrix),
           norm_meta_alpha(Metacommunity(pop, Z), qs))
@deprecate(subcommunityalpha(pop::AbstractArray, qs, Z::AbstractMatrix),
           raw_meta_alpha(Metacommunity(pop, Z), qs))
@deprecate(subcommunitybetabar(pop::AbstractArray, qs, Z::AbstractMatrix),
           norm_meta_beta(Metacommunity(pop, Z), qs))
@deprecate(subcommunitybeta(pop::AbstractArray, qs, Z::AbstractMatrix),
           raw_meta_beta(Metacommunity(pop, Z), qs))
@deprecate(subcommunityrhobar(pop::AbstractArray, qs, Z::AbstractMatrix),
           norm_meta_rho(Metacommunity(pop, Z), qs))
@deprecate(subcommunityrho(pop::AbstractArray, qs, Z::AbstractMatrix),
           raw_meta_rho(Metacommunity(pop, Z), qs))
@deprecate(subcommunitygammabar(pop::AbstractArray, qs, Z::AbstractMatrix),
           meta_gamma(Metacommunity(pop, Z), qs))
@deprecate(subcommunitygamma(pop::AbstractArray, qs, Z::AbstractMatrix),
           meta_gamma(Metacommunity(pop, Z), qs))

@deprecate(subcommunityalphabar(pop::AbstractArray, qs),
           norm_meta_alpha(Metacommunity(pop), qs))
@deprecate(subcommunityalpha(pop::AbstractArray, qs),
           raw_meta_alpha(Metacommunity(pop), qs))
@deprecate(subcommunitybetabar(pop::AbstractArray, qs),
           norm_meta_beta(Metacommunity(pop), qs))
@deprecate(subcommunitybeta(pop::AbstractArray, qs),
           raw_meta_beta(Metacommunity(pop), qs))
@deprecate(subcommunityrhobar(pop::AbstractArray, qs),
           norm_meta_rho(Metacommunity(pop), qs))
@deprecate(subcommunityrho(pop::AbstractArray, qs),
           raw_meta_rho(Metacommunity(pop), qs))
@deprecate(subcommunitygammabar(pop::AbstractArray, qs),
           meta_gamma(Metacommunity(pop), qs))
@deprecate(subcommunitygamma(pop::AbstractArray, qs),
           meta_gamma(Metacommunity(pop), qs))

## Deprecate all-in-one names, as we divide calculation into type of
## diversity and scale
@deprecate(normsubalpha(pop::AbstractArray, Z, qs),
           norm_sub_alpha(Metacommunity(pop, Z), qs))
@deprecate(rawsubalpha(pop::AbstractArray, Z, qs),
           raw_sub_alpha(Metacommunity(pop, Z), qs))
@deprecate(normsubbeta(pop::AbstractArray, Z, qs),
           norm_sub_beta(Metacommunity(pop, Z), qs))
@deprecate(rawsubbeta(pop::AbstractArray, Z, qs),
           raw_sub_beta(Metacommunity(pop, Z), qs))
@deprecate(normsubrho(pop::AbstractArray, Z, qs),
           norm_sub_rho(Metacommunity(pop, Z), qs))
@deprecate(rawsubrho(pop::AbstractArray, Z, qs),
           raw_sub_rho(Metacommunity(pop, Z), qs))
@deprecate(subgamma(pop::AbstractArray, Z, qs),
           sub_gamma(Metacommunity(pop, Z), qs))
@deprecate(normmetaalpha(pop::AbstractArray, Z, qs),
           norm_meta_alpha(Metacommunity(pop, Z), qs))
@deprecate(rawmetaalpha(pop::AbstractArray, Z, qs),
           raw_meta_alpha(Metacommunity(pop, Z), qs))
@deprecate(normmetabeta(pop::AbstractArray, Z, qs),
           norm_meta_beta(Metacommunity(pop, Z), qs))
@deprecate(rawmetabeta(pop::AbstractArray, Z, qs),
           raw_meta_beta(Metacommunity(pop, Z), qs))
@deprecate(normmetarho(pop::AbstractArray, Z, qs),
           norm_meta_rho(Metacommunity(pop, Z), qs))
@deprecate(rawmetarho(pop::AbstractArray, Z, qs),
           raw_meta_rho(Metacommunity(pop, Z), qs))
@deprecate(metagamma(pop::AbstractArray, Z, qs),
           meta_gamma(Metacommunity(pop, Z), qs))

@deprecate(normsubalpha(pop::AbstractArray, qs),
           norm_sub_alpha(Metacommunity(pop), qs))
@deprecate(rawsubalpha(pop::AbstractArray, qs),
           raw_sub_alpha(Metacommunity(pop), qs))
@deprecate(normsubbeta(pop::AbstractArray, qs),
           norm_sub_beta(Metacommunity(pop), qs))
@deprecate(rawsubbeta(pop::AbstractArray, qs),
           raw_sub_beta(Metacommunity(pop), qs))
@deprecate(normsubrho(pop::AbstractArray, qs),
           norm_sub_rho(Metacommunity(pop), qs))
@deprecate(rawsubrho(pop::AbstractArray, qs),
           raw_sub_rho(Metacommunity(pop), qs))
@deprecate(subgamma(pop::AbstractArray, qs),
           sub_gamma(Metacommunity(pop), qs))
@deprecate(normmetaalpha(pop::AbstractArray, qs),
           norm_meta_alpha(Metacommunity(pop), qs))
@deprecate(rawmetaalpha(pop::AbstractArray, qs),
           raw_meta_alpha(Metacommunity(pop), qs))
@deprecate(normmetabeta(pop::AbstractArray, qs),
           norm_meta_beta(Metacommunity(pop), qs))
@deprecate(rawmetabeta(pop::AbstractArray, qs),
           raw_meta_beta(Metacommunity(pop), qs))
@deprecate(normmetarho(pop::AbstractArray, qs),
           norm_meta_rho(Metacommunity(pop), qs))
@deprecate(rawmetarho(pop::AbstractArray, qs),
           raw_meta_rho(Metacommunity(pop), qs))
@deprecate(metagamma(pop::AbstractArray, qs),
           meta_gamma(Metacommunity(pop), qs))

## Deprecate anything related to ϵ as it has been replaced by ρ̄
@deprecate(subcommunityepsilon, Dρ̄)
@deprecate(Dϵ, Dρ̄)
@deprecate(ecosystemE, ecosystemRbar)
@deprecate(supercommunityE, supercommunityRbar)
@deprecate(metacommunityE, metacommunityRbar)
@deprecate(DE, metacommunityE)

"""
The **Diversity.Ecology** module replicates old ecological
diversity measures and generalised versions of them that relate to our
general measures of alpha, beta and gamma diversity at subcommunity
and metacommunity levels. The generalisations of the richness, Shannon
and Simpson are the only standard measures we are aware of whose
subcommunity components sum directly to the corresponding ecosystem
measure (although note that Simpson's index decreases for increased
diversity, so small components are more diverse).
"""
module Ecology

include("Ecology.jl")
export generalisedrichness, richness
export generalisedshannon, shannon
export generalisedsimpson, simpson
export generalisedjaccard, jaccard

end # sub-module Ecology

"""
Lou Jost's
[diversity](http://dx.doi.org/10.1111/j.2006.0030-1299.14714.x)
[measures](http://www.esajournals.org/doi/abs/10.1890/06-1736.1) are
found in the **Diversity.Jost** module.
"""
module Jost

include("Jost.jl")
export jostbeta, jostalpha
@deprecate(jostα, jostalpha)
@deprecate(jostβ, jostbeta)

end # sub-module Jost

"""
[Hill numbers](http://www.jstor.org/stable/1934352) are found in the
**Diversity.Hill** package.
"""
module Hill

include("Hill.jl")
export hillnumber

end # sub-module Hill

end # module
