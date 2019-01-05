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
# Base types
export AbstractPartition, AbstractTypes, AbstractMetacommunity

# Base class and functions required for each partition
export _getsubcommunitynames # required
export _countsubcommunities  # optional

# Base class and functions required for each type
export _gettypenames, _calcsimilarity, _getscale      # required
export _counttypes, _calcabundance, _calcordinariness # optional
export _getdiversityname, _addedoutputcols, _getaddedoutput # optional

# Base class and functions required for each metacommunity
export _gettypes, _getpartition, _getabundance  # required
export _getmetaabundance, _getweight            # optional
export _getordinariness!, _getmetaordinariness! # optional
export _hassimilarity

# Function with minimal default implementation for types and partitions
export floattypes
# Functions with standard implementations
export typematch, mcmatch
end

# Add in our user interface and the EcoBase higher level interface
include("Interface.jl")
export gettypes, getpartition
export counttypes, countsubcommunities
export gettypenames, getsubcommunitynames
export getabundance, getmetaabundance, getweight
export getordinariness!, getmetaordinariness!
export hassimilarity, calcsimilarity # Needed because it sometimes doesn't exist unless requested
export getdiversityname, addedoutputcols, getaddedoutput

# Inheritance from EcoBase
include("EcoBase.jl")

include("Iterators.jl")
export TypeIterator, SubcommunityIterator

include("Types.jl")
export UniqueTypes, Species, Taxonomy, GeneralTypes

include("Partition.jl")
export Subcommunities, Onecommunity

include("Metacommunity.jl")
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

# Does PhyloTypes need to exist already?
using Requires
@static if VERSION < v"0.7.0-"
    @require Phylo begin
        println("Creating Diversity to Phylo interface...")
        include("Phylogenetics.jl")
        export PhyloTypes
    end
else
    function __init__()
        @require Phylo="aea672f4-3940-5932-aa44-993d1c3ff149" begin
            println("Creating Diversity to Phylo interface...")
            include("Phylogenetics.jl")
            export PhyloTypes
        end
    end
end

include("GeneralisedDiversities.jl")
export diversity
export norm_sub_alpha, raw_sub_alpha, norm_sub_beta, raw_sub_beta
export norm_sub_rho, raw_sub_rho, sub_gamma
export norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta
export norm_meta_rho, raw_meta_rho, meta_gamma

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

# Path into package
path(path...; dir::String = "test") = joinpath(@__DIR__, "..", dir, path...)

end # module
