module Diversity

export powermean
export qD, qDZ

export ᾱ, communityalphabar, α, communityalpha, A, ecosystemA, Ā, ecosystemAbar
export β̄, communitybetabar, β, communitybeta, B, ecosystemB, B̄, ecosystemBbar
export γ̄, communitygammabar, γ, communitygamma, G, ecosystemG, Ḡ, ecosystemGbar

export generalisedrichness, richness
export generalisedshannon, shannon
export generalisedsimpson, simpson

include("EffectiveNumbers.jl")
include("GeneralisedDiversities.jl")
#include("CommunityContributions.jl")
include("HistoricalDiversities.jl")

end # module
