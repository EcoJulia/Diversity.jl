module Diversity
using Docile
@docstrings [ :manual => ["../doc/diversity.md"] ]

include("EffectiveNumbers.jl")
export qD, qDZ

include("GeneralisedDiversities.jl")
export ᾱ, subcommunityalphabar, α, subcommunityalpha, A, ecosystemA, Ā, ecosystemAbar
export β̄, subcommunitybetabar, β, subcommunitybeta, B, ecosystemB, B̄, ecosystemBbar
export ρ̄, subcommunityrhobar, ρ, subcommunityrho, R, ecosystemR, R̄, ecosystemRbar
export ϵ, subcommunityevenness, subcommunityredundancy, ecosystemevenness, ecosystemredundancy
export γ̄, subcommunitygammabar, γ, subcommunitygamma, G, ecosystemG, Ḡ, ecosystemGbar
export diversity

include("CommunityContributions.jl")
export contributions

module Ecology
using Docile
@docstrings [ :manual => ["../doc/ecology.md"] ]

include("Ecology.jl")
export generalisedrichness, richness
export generalisedshannon, shannon
export generalisedsimpson, simpson
export generalisedjaccard, jaccard

end # sub-module Ecology

module Jost
using Docile
@docstrings [ :manual => ["../doc/jost.md"] ]

include("Jost.jl")
export jostD, jostbeta, jostβ

end # sub-module Jost

module Hill
using Docile
@docstrings [ :manual => ["../doc/hill.md"] ]

include("Hill.jl")
export hillnumber

end # sub-module Hill

end # module
