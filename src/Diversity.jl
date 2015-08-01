module Diversity
if VERSION < v"0.4-"
    using Docile
end

if VERSION < v"0.4-"
   @docstrings(manual = ["../doc/diversity.md"])
end

#include("Collection.jl")
#export Collection, Ecosystem, Community
#export Subcommunity, Onecommunity
#export Unique, Species, Taxonomy

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
if VERSION < v"0.4-"
   @docstrings(manual = ["../doc/ecology.md"])
end

include("Ecology.jl")
export generalisedrichness, richness
export generalisedshannon, shannon
export generalisedsimpson, simpson
export generalisedjaccard, jaccard

end # sub-module Ecology

module Jost
using Docile
if VERSION < v"0.4-"
   @docstrings(manual = ["../doc/jost.md"])
end

include("Jost.jl")
export jostbeta, jostβ, jostalpha, jostα

end # sub-module Jost

module Hill
using Docile
if VERSION < v"0.4-"
   @docstrings(manual = ["../doc/hill.md"])
end

include("Hill.jl")
export hillnumber

end # sub-module Hill

## Make sure that Lexicon is loading so inline REPL documentation works
using Lexicon

end # module
