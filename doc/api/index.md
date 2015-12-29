# API-INDEX


## MODULE: Diversity

---

## Functions [Exported]

[Diversity.Community](Diversity.md#function__community.1)  ### Community type, representing a single community

[Diversity.Ecosystem](Diversity.md#function__ecosystem.1)  ### Ecosystem type, representing an ecosystem of multiple subcommunities

---

## Methods [Exported]

[diversity{S<:AbstractFloat, T<:Diversity.Similarity}(kind::Tuple{Function, Set{Symbol}},  proportions::Array{S<:AbstractFloat, 2},  qs,  sim::T<:Diversity.Similarity)](Diversity.md#method__diversity.1)  ### Calculates subcommunity and supercommunity diversities

[diversity{S<:AbstractFloat}(kind::Tuple{Function, Set{Symbol}},  proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__diversity.2)  ### Calculates subcommunity and supercommunity diversities

[diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__diversity.3)  ### Calculates subcommunity and supercommunity diversities

[diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  sim,  returnsupercommunity::Bool)](Diversity.md#method__diversity.4)  ### Calculates subcommunity and supercommunity diversities

[diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  sim,  returnsupercommunity::Bool,  returnsubcommunity::Bool)](Diversity.md#method__diversity.5)  ### Calculates subcommunity and supercommunity diversities

[diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  sim,  returnsupercommunity::Bool,  returnsubcommunity::Bool,  returnweights::Bool)](Diversity.md#method__diversity.6)  ### Calculates subcommunity and supercommunity diversities

[qDZ{S<:AbstractFloat, T<:Diversity.Similarity}(proportions::Array{S<:AbstractFloat, 1},  qs,  sim::T<:Diversity.Similarity)](Diversity.md#method__qdz.1)  ### Calculates Leinster-Cobbold / similarity-sensitive diversity

[qDZ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 1},  qs)](Diversity.md#method__qdz.2)  ### Calculates Leinster-Cobbold / similarity-sensitive diversity

[qD{S<:Real}(proportions::Array{S<:Real, 1},  qs)](Diversity.md#method__qd.1)  ### Calculates Hill / naive-similarity diversity

[subcommunityalphabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__subcommunityalphabar.1)  ### Normalised similarity-sensitive subcommunity alpha diversity)

[subcommunityalphabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__subcommunityalphabar.2)  ### Normalised similarity-sensitive subcommunity alpha diversity)

[subcommunityalpha{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__subcommunityalpha.1)  ### Raw similarity-sensitive subcommunity alpha diversity / naive-community diversity

[subcommunityalpha{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__subcommunityalpha.2)  ### Raw similarity-sensitive subcommunity alpha diversity / naive-community diversity

[subcommunitybetabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__subcommunitybetabar.1)  ### Normalised similarity-sensitive subcommunity beta diversity

[subcommunitybetabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__subcommunitybetabar.2)  ### Normalised similarity-sensitive subcommunity beta diversity

[subcommunitybeta{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__subcommunitybeta.1)  ### Raw similarity-sensitive subcommunity beta diversity / distinctiveness / concentration

[subcommunitybeta{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__subcommunitybeta.2)  ### Raw similarity-sensitive subcommunity beta diversity / distinctiveness / concentration

[subcommunitygammabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__subcommunitygammabar.1)  ### Normalised similarity-sensitive subcommunity gamma diversity

[subcommunitygammabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__subcommunitygammabar.2)  ### Normalised similarity-sensitive subcommunity gamma diversity

[subcommunitygamma{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__subcommunitygamma.1)  ### Raw similarity-sensitive subcommunity gamma diversity

[subcommunitygamma{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__subcommunitygamma.2)  ### Raw similarity-sensitive subcommunity gamma diversity

[subcommunityrhobar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__subcommunityrhobar.1)  ### Normalised similarity-sensitive subcommunity representativeness

[subcommunityrhobar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__subcommunityrhobar.2)  ### Normalised similarity-sensitive subcommunity representativeness

[subcommunityrho{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__subcommunityrho.1)  ### Raw similarity-sensitive subcommunity redundancy

[subcommunityrho{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__subcommunityrho.2)  ### Raw similarity-sensitive subcommunity redundancy

[supercommunityAbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__supercommunityabar.1)  ### Normalised similarity-sensitive supercommunity alpha diversity

[supercommunityAbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__supercommunityabar.2)  ### Normalised similarity-sensitive supercommunity alpha diversity

[supercommunityA{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__supercommunitya.1)  ### Raw similarity-sensitive supercommunity alpha diversity / naive-community diversity

[supercommunityA{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__supercommunitya.2)  ### Raw similarity-sensitive supercommunity alpha diversity / naive-community diversity

[supercommunityBbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__supercommunitybbar.1)  ### Normalised similarity-sensitive supercommunity beta diversity / effective number of communities

[supercommunityBbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__supercommunitybbar.2)  ### Normalised similarity-sensitive supercommunity beta diversity / effective number of communities

[supercommunityB{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__supercommunityb.1)  ### Raw similarity-sensitive supercommunity beta diversity / distinctiveness / concentration

[supercommunityB{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__supercommunityb.2)  ### Raw similarity-sensitive supercommunity beta diversity / distinctiveness / concentration

[supercommunityGbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__supercommunitygbar.1)  ### Normalised similarity-sensitive supercommunity gamma diversity

[supercommunityGbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__supercommunitygbar.2)  ### Normalised similarity-sensitive supercommunity gamma diversity

[supercommunityG{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__supercommunityg.1)  ### Raw similarity-sensitive supercommunity gamma diversity

[supercommunityG{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__supercommunityg.2)  ### Raw similarity-sensitive supercommunity gamma diversity

[supercommunityRbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__supercommunityrbar.1)  ### Normalised similarity-sensitive supercommunity representativeness

[supercommunityRbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__supercommunityrbar.2)  ### Normalised similarity-sensitive supercommunity representativeness

[supercommunityR{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__supercommunityr.1)  ### Raw similarity-sensitive supercommunity redundancy

[supercommunityR{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim)](Diversity.md#method__supercommunityr.2)  ### Raw similarity-sensitive supercommunity redundancy

---

## Types [Exported]

[Diversity.Collection{S<:Diversity.Similarity, P<:Diversity.Partition, FP<:AbstractFloat}](Diversity.md#type__collection.1)  ### Collection type, representing a collection of one or more subcommunities

[Diversity.GeneralSimilarity{S<:AbstractFloat}](Diversity.md#type__generalsimilarity.1)  ### A general matrix-based Similarity subtype

[Diversity.Onecommunity](Diversity.md#type__onecommunity.1)  ### Partition type allowing only one subcommunity

[Diversity.Subcommunity](Diversity.md#type__subcommunity.1)  ### Partition type with multiple subccomunities

[Diversity.Taxonomy](Diversity.md#type__taxonomy.1)  ### A subtype of Similarity with similarity between related taxa

[Diversity.Unique](Diversity.md#type__unique.1)  ### A subtype of Similarity where all individuals are completely distinct

---

## Typealiass [Exported]

[Species](Diversity.md#typealias__species.1)  ### A subtype of Similarity where all species are completely distinct

---

## Methods [Internal]

[call{S<:AbstractFloat}(::Type{Diversity.GeneralSimilarity{S<:AbstractFloat}},  z::Array{S<:AbstractFloat, 2})](Diversity.md#method__call.1)  ### Constructor for GeneralSimilarity

[contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.md#method__contributions.1)  ### Calculate diversity contributions from subcommunities

[contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  perindividual::Bool)](Diversity.md#method__contributions.2)  ### Calculate diversity contributions from subcommunities

[contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  perindividual::Bool,  Z::Array{S<:AbstractFloat, 2})](Diversity.md#method__contributions.3)  ### Calculate diversity contributions from subcommunities

[contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  perindividual::Bool,  Z::Array{S<:AbstractFloat, 2},  returnsupercommunity::Bool)](Diversity.md#method__contributions.4)  ### Calculate diversity contributions from subcommunities

[contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  perindividual::Bool,  Z::Array{S<:AbstractFloat, 2},  returnsupercommunity::Bool,  returnsubcommunity::Bool)](Diversity.md#method__contributions.5)  ### Calculate diversity contributions from subcommunities

[contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  perindividual::Bool,  Z::Array{S<:AbstractFloat, 2},  returnsupercommunity::Bool,  returnsubcommunity::Bool,  returnweights::Bool)](Diversity.md#method__contributions.6)  ### Calculate diversity contributions from subcommunities

[powermean{S<:AbstractFloat}(values::Array{S<:AbstractFloat, 1})](Diversity.md#method__powermean.1)  ### Calculates the weighted powermean of a series of numbers

[powermean{S<:AbstractFloat}(values::Array{S<:AbstractFloat, 1},  order::S<:AbstractFloat)](Diversity.md#method__powermean.2)  ### Calculates the weighted powermean of a series of numbers

[powermean{S<:AbstractFloat}(values::Array{S<:AbstractFloat, 1},  order::S<:AbstractFloat,  weights::Array{S<:AbstractFloat, 1})](Diversity.md#method__powermean.3)  ### Calculates the weighted powermean of a series of numbers

---

## Types [Internal]

[Diversity.Partition](Diversity.md#type__partition.1)  ### Abstract Partition supertype for all partitioning types

[Diversity.Similarity](Diversity.md#type__similarity.1)  ### Abstract Similarity supertype for all similarity measures

## MODULE: Diversity.Jost

---

## Methods [Exported]

[jostalpha{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.Jost.md#method__jostalpha.1)  ### Calculates Jost's alpha diversity

[jostbeta{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs)](Diversity.Jost.md#method__jostbeta.1)  ### Calculates Jost's beta diversity

## MODULE: Diversity.Ecology

---

## Methods [Exported]

[generalisedjaccard(proportions::Array{T, 2},  qs)](Diversity.Ecology.md#method__generalisedjaccard.1)  ### Calculate a generalised version of the Jaccard index

[generalisedjaccard(proportions::Array{T, 2},  qs,  Z::Array{T, 2})](Diversity.Ecology.md#method__generalisedjaccard.2)  ### Calculate a generalised version of the Jaccard index

[generalisedrichness{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2})](Diversity.Ecology.md#method__generalisedrichness.1)  ### Calculate a generalised version of richness

[generalisedrichness{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  Z::Array{S<:AbstractFloat, 2})](Diversity.Ecology.md#method__generalisedrichness.2)  ### Calculate a generalised version of richness

[generalisedshannon{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2})](Diversity.Ecology.md#method__generalisedshannon.1)  ### Calculate a generalised version of Shannon entropy

[generalisedshannon{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  Z::Array{S<:AbstractFloat, 2})](Diversity.Ecology.md#method__generalisedshannon.2)  ### Calculate a generalised version of Shannon entropy

[generalisedsimpson{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2})](Diversity.Ecology.md#method__generalisedsimpson.1)  ### Calculate a generalised version of Simpson's index

[generalisedsimpson{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  Z::Array{S<:AbstractFloat, 2})](Diversity.Ecology.md#method__generalisedsimpson.2)  ### Calculate a generalised version of Simpson's index

[jaccard{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2})](Diversity.Ecology.md#method__jaccard.1)  ### Calculate the Jaccard index

[richness{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2})](Diversity.Ecology.md#method__richness.1)  ### Calculate species richness of populations

[shannon{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2})](Diversity.Ecology.md#method__shannon.1)  ### Calculate Shannon entropy of populations

[simpson{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2})](Diversity.Ecology.md#method__simpson.1)  ### Calculate Simpson's index

## MODULE: Diversity.Hill

---

## Methods [Exported]

[hillnumber(proportions,  qs)](Diversity.Hill.md#method__hillnumber.1)  ### Calculates Hill numbers

