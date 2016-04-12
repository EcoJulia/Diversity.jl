__precompile__()

"$((VERSION < v"0.5-"? readall : readstring)(joinpath(dirname(@__FILE__), "../doc/diversity.md")))"
module Diversity
using Compat

include("Supercommunity.jl")
export Partition, Subcommunities, Onecommunity
export Similarity, Unique, Species, Taxonomy, MatrixSimilarity
export Supercommunity, Ecosystem, Community
export getAbundances, getWeights
export getSimilarityMatrix, getOrdinariness!, getSuperOrdinariness!

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

"$(@compat readstring(joinpath(dirname(@__FILE__), "../doc/ecology.md")))"
module Ecology
using Compat

include("Ecology.jl")
export generalisedrichness, richness
export generalisedshannon, shannon
export generalisedsimpson, simpson
export generalisedjaccard, jaccard

end # sub-module Ecology

"$(@compat readstring(joinpath(dirname(@__FILE__), "../doc/jost.md")))"
module Jost
using Compat

include("Jost.jl")
export jostbeta, jostβ, jostalpha, jostα

end # sub-module Jost

"$(@compat readstring(joinpath(dirname(@__FILE__), "../doc/hill.md")))"
module Hill
using Compat

include("Hill.jl")
export hillnumber

end # sub-module Hill

end # module
