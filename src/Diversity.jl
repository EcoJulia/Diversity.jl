__precompile__()

"""
The main **Diversity** module provides basic numbers-equivalent
diversity measures (described in
[Hill, 1973](http://www.jstor.org/stable/1934352)),
similarity-sensitive diversity measures (generalised from Hill, and
described in
[Leinster and Cobbold, 2012](http://www.esajournals.org/doi/abs/10.1890/10-2402.1)),
and related alpha, beta and gamma diversity measures at the level of
the supercommunity and its component subcommunities (generalised in
turn from Leinster and Cobbold, and described in
[Reeve et al, 2014](http://arxiv.org/abs/1404.6520)). The diversity
functions exist both with unicode names (e.g. ```ᾱ()```), which are
not automatically exported (as we feel they are too short) and with
matching longer ASCII names (e.g. ```NormalisedAlpha()```), which are.
We also provide functions to calculate appropriate
```subcommunityDiversity()``` and ```supercommunityDiversity()```
values for each measure, a general ```diversity()``` function for
extract any diversity measure at a series of scales.
"""
module Diversity

include("Supercommunity.jl")
export Subcommunities, Onecommunity
export Similarity, Unique, Species, Taxonomy, MatrixSimilarity
export Supercommunity, Ecosystem, SingleCommunity
export getabundance, getweight
export getsimilarity, getordinariness!, getsuperordinariness!
export inddiv, subdiv, superdiv

include("EffectiveNumbers.jl")
export qD, qDZ

include("DiversityMeasure.jl")
export DiversityLevel
export individualDiversity, subcommunityDiversity, supercommunityDiversity

export RawAlpha, NormalisedAlpha
export RawBeta, NormalisedBeta, RawRho, NormalisedRho
export Gamma
## We do not export ᾱ, α, β̄, β, ρ̄, ρ, γ̄, γ as they're too short

export getName, getASCIIName, getFullName

include("GeneralisedDiversities.jl")
export diversity

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
@deprecate(ecosystemA, supercommunityA)
@deprecate(ecosystemAbar, supercommunityAbar)
@deprecate(ecosystemB, supercommunityB)
@deprecate(ecosystemBbar, supercommunityBbar)
@deprecate(ecosystemR, supercommunityR)
@deprecate(ecosystemRbar, supercommunityRbar)
@deprecate(ecosystemG, supercommunityG)
@deprecate(ecosystemGbar, supercommunityGbar)

## Deprecate all-in-one names, as we divide calculation into type of
## diversity and scale
@deprecate(subcommunityalphabar(pop, Z, qs),
           subcommunityDiversity(ᾱ(Ecosystem(pop, Z)), qs))
@deprecate(subcommunityalpha(pop, Z, qs),
           subcommunityDiversity(α(Ecosystem(pop, Z)), qs))
@deprecate(subcommunitybetabar(pop, Z, qs),
           subcommunityDiversity(β̄(Ecosystem(pop, Z)), qs))
@deprecate(subcommunitybeta(pop, Z, qs),
           subcommunityDiversity(β(Ecosystem(pop, Z)), qs))
@deprecate(subcommunityrhobar(pop, Z, qs),
           subcommunityDiversity(ρ̄(Ecosystem(pop, Z)), qs))
@deprecate(subcommunityrho(pop, Z, qs),
           subcommunityDiversity(ρ(Ecosystem(pop, Z)), qs))
@deprecate(subcommunitygammabar(pop, Z, qs),
           subcommunityDiversity(γ̄(Ecosystem(pop, Z)), qs))
@deprecate(subcommunitygamma(pop, Z, qs),
           subcommunityDiversity(γ(Ecosystem(pop, Z)), qs))
@deprecate(supercommunityAbar(pop, Z, qs),
           supercommunityDiversity(ᾱ(Ecosystem(pop, Z)), qs))
@deprecate(supercommunityA(pop, Z, qs),
           supercommunityDiversity(α(Ecosystem(pop, Z)), qs))
@deprecate(supercommunityBbar(pop, Z, qs),
           supercommunityDiversity(β̄(Ecosystem(pop, Z)), qs))
@deprecate(supercommunityB(pop, Z, qs),
           supercommunityDiversity(β(Ecosystem(pop, Z)), qs))
@deprecate(supercommunityRbar(pop, Z, qs),
           supercommunityDiversity(ρ̄(Ecosystem(pop, Z)), qs))
@deprecate(supercommunityR(pop, Z, qs),
           supercommunityDiversity(ρ(Ecosystem(pop, Z)), qs))
@deprecate(supercommunityG(pop, Z, qs),
           supercommunityDiversity(γ(Ecosystem(pop, Z)), qs))

@deprecate(subcommunityalphabar(pop, qs),
           subcommunityDiversity(ᾱ(Ecosystem(pop)), qs))
@deprecate(subcommunityalpha(pop, qs),
           subcommunityDiversity(α(Ecosystem(pop)), qs))
@deprecate(subcommunitybetabar(pop, qs),
           subcommunityDiversity(β̄(Ecosystem(pop)), qs))
@deprecate(subcommunitybeta(pop, qs),
           subcommunityDiversity(β(Ecosystem(pop)), qs))
@deprecate(subcommunityrhobar(pop, qs),
           subcommunityDiversity(ρ̄(Ecosystem(pop)), qs))
@deprecate(subcommunityrho(pop, qs),
           subcommunityDiversity(ρ(Ecosystem(pop)), qs))
@deprecate(subcommunitygammabar(pop, qs),
           subcommunityDiversity(γ̄(Ecosystem(pop)), qs))
@deprecate(subcommunitygamma(pop, qs),
           subcommunityDiversity(γ(Ecosystem(pop)), qs))
@deprecate(supercommunityAbar(pop, qs),
           supercommunityDiversity(ᾱ(Ecosystem(pop)), qs))
@deprecate(supercommunityA(pop, qs),
           supercommunityDiversity(α(Ecosystem(pop)), qs))
@deprecate(supercommunityBbar(pop, qs),
           supercommunityDiversity(β̄(Ecosystem(pop)), qs))
@deprecate(supercommunityB(pop, qs),
           supercommunityDiversity(β(Ecosystem(pop)), qs))
@deprecate(supercommunityRbar(pop, qs),
           supercommunityDiversity(ρ̄(Ecosystem(pop)), qs))
@deprecate(supercommunityR(pop, qs),
           supercommunityDiversity(ρ(Ecosystem(pop)), qs))
@deprecate(supercommunityG(pop, qs),
           supercommunityDiversity(γ(Ecosystem(pop)), qs))

## Deprecate anything related to ϵ as it has been replaced by ρ̄
@deprecate(subcommunityepsilon, Dρ̄)
@deprecate(Dϵ, Dρ̄)
@deprecate(ecosystemE, supercommunityRbar)
@deprecate(supercommunityE, supercommunityRbar)
@deprecate(DE, supercommunityRbar)

"""
The **Diversity.Ecology** module replicates old ecological
diversity measures and generalised versions of them that relate to our
general measures of alpha, beta and gamma diversity at subcommunity
and supercommunity levels. The generalisations of the richness, Shannon
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
export jostbeta, jostβ, jostalpha, jostα

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
