module Diversity

include("EffectiveNumbers.jl")
export powermean
export qD, qDZ

include("GeneralisedDiversities.jl")
export ᾱ, communityalphabar, α, communityalpha, A, ecosystemA, Ā, ecosystemAbar
export β̄, communitybetabar, β, communitybeta, B, ecosystemB, B̄, ecosystemBbar
export γ̄, communitygammabar, γ, communitygamma, G, ecosystemG, Ḡ, ecosystemGbar

include("CommunityContributions.jl")
export contributions

include("HistoricalDiversities.jl")
export generalisedrichness, richness
export generalisedshannon, shannon
export generalisedsimpson, simpson

end # module
