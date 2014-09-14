module Diversity

include("EffectiveNumbers.jl")
export qD, qDZ

include("GeneralisedDiversities.jl")
export ᾱ, subcommunityalphabar, α, subcommunityalpha, A, ecosystemA, Ā, ecosystemAbar
export β̄, subcommunitybetabar, β, subcommunitybeta, B, ecosystemB, B̄, ecosystemBbar
export γ̄, subcommunitygammabar, γ, subcommunitygamma, G, ecosystemG, Ḡ, ecosystemGbar
export diversity

include("CommunityContributions.jl")
export contributions

module Compatibility

include("HistoricalDiversities.jl")
export generalisedrichness, richness
export generalisedshannon, shannon
export generalisedsimpson, simpson
export generalisedjaccard, jaccard

end # sub-module Compatibility

module Jost

include("Jost.jl")
export jostD, jostbeta, jostβ

end # sub-module Jost

module Hill

include("Hill.jl")
export hillnumber

end # sub-module Hill

end # module
