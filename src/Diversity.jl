module Diversity
using Docile
@docstrings(manual = ["../doc/diversity.md"])

include("Collection.jl")
export Collection, Ecosystem, Community
export Subcommunity, Onecommunity
export Unique, Species, Taxonomy

include("EffectiveNumbers.jl")
export qD, qDZ

include("GeneralisedDiversities.jl")
export subcommunityalphabar, subcommunityalpha, ecosystemA, ecosystemAbar
export subcommunitybetabar, subcommunitybeta, ecosystemB, ecosystemBbar
export subcommunitygammabar, subcommunitygamma, ecosystemG, ecosystemGbar
export Dᾱ, Dα, Dβ̄, Dβ, Dϵ, Dρ̄, Dρ, Dγ̄, Dγ
export DĀ, DA, DB̄, DB, DE, DR̄, DR, DḠ, DG
export diversity

include("CommunityContributions.jl")
export contributions

module Ecology
using Docile
@docstrings(manual = ["../doc/ecology.md"])

include("Ecology.jl")
export generalisedrichness, richness
export generalisedshannon, shannon
export generalisedsimpson, simpson
export generalisedjaccard, jaccard

end # sub-module Ecology

module Jost
using Docile
@docstrings(manual = ["../doc/jost.md"])

include("Jost.jl")
export jostD, jostbeta, jostβ

end # sub-module Jost

module Hill
using Docile
@docstrings(manual = ["../doc/hill.md"])

include("Hill.jl")
export hillnumber

end # sub-module Hill

end # module
