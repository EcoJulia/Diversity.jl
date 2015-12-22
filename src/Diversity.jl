module Diversity
using Docile

"""
!!summary(Package for measuring diversity for the Julia Language)

!!include(../doc/diversity.md)
"""
Diversity

include("Collection.jl")
export Collection, Ecosystem, Community
export Subcommunity, Onecommunity
export Unique, Species, Taxonomy

include("EffectiveNumbers.jl")
export qD, qDZ

include("GeneralisedDiversities.jl")
export subcommunityalphabar, subcommunityalpha, supercommunityA, supercommunityAbar, ecosystemA, ecosystemAbar
export subcommunitybetabar, subcommunitybeta, supercommunityB, supercommunityBbar, ecosystemB, ecosystemBbar
export subcommunityrhobar, subcommunityrho, supercommunityR, supercommunityRbar, ecosystemR, ecosystemRbar
export subcommunitygammabar, subcommunitygamma, supercommunityG, supercommunityGbar, ecosystemG, ecosystemGbar
export Dᾱ, Dα, Dβ̄, Dβ, Dϵ, Dρ̄, Dρ, Dγ̄, Dγ
export DĀ, DA, DB̄, DB, DE, DR̄, DR, DḠ, DG
export diversity

include("CommunityContributions.jl")
export contributions

module Ecology
using Docile

"""
!!summary(Sub-package of Diversity for ecological diversity measures)

!!include(../doc/ecology.md)
"""
Ecology

include("Ecology.jl")
export generalisedrichness, richness
export generalisedshannon, shannon
export generalisedsimpson, simpson
export generalisedjaccard, jaccard

end # sub-module Ecology

module Jost
using Docile

"""
!!summary(Sub-package of Diversity for Jost's diversity measures)

!!include(../doc/jost.md)
"""
Jost

include("Jost.jl")
export jostbeta, jostβ, jostalpha, jostα

end # sub-module Jost

module Hill
using Docile

"""
!!summary(Sub-package of Diversity for Hill numbers)

!!include(../doc/hill.md)
"""
Hill

include("Hill.jl")
export hillnumber

end # sub-module Hill

## Make sure that Lexicon is loading so inline REPL documentation works
using Lexicon

end # module
