# Diversity.Ecology

In the **Diversity.Ecology** submodule, we replicate old ecological
diversity measures and generalised versions of them that relate to our
general measures of alpha, beta and gamma diversity at subcommunity
and ecosystem measures. The generalisations of the richness, Shannon
and Simpson are the only standard measures we are aware of whose
subcommunity components sum directly to the corresponding ecosystem
measure (although note that Simpson's index decreases for increased
diversity, so small components are more diverse).

## Usage

Accessing the functionality in the package is simple. Note that the submodule
provides the ecological measures themselves, while `Metacommunity` and the
[`DiversityLevel`](@ref)s come from `Diversity`, so both are loaded here:

```@repl ecology
using Diversity
using Diversity.Ecology
community = [10, 20, 20]
community = community ./ sum(community)
simpson(community)
shannon(community)
richness(community)
```

Two subcommunities can be compared with the Jaccard index, either directly or -
since it is a special case of our general measures - through a `Metacommunity`,
with or without a similarity matrix:

```@repl ecology
using LinearAlgebra
ecosystem = [2 2 0; 0 2 2]'
ecosystem = ecosystem ./ sum(ecosystem)
jaccard(ecosystem)
generalisedjaccard(Metacommunity(ecosystem))
generalisedjaccard(ecosystem, Matrix(1.0I, 3, 3))
```

Pielou's evenness measures how equally the individuals are spread across the
types, from zero to one:

```@repl ecology
pielou([0.7, 0.2, 0.1])
communitymat = [10 20 30 20 0;   # 5 subcommunities (columns), 6 species (rows)
                10  0 50 80 10;
                60 10 90  0  0;
                10 10 10 10 10;
                70 70 70 70 70;
                10  0  0 90  0]
Diversity.Ecology.generalisedpielou(subcommunityDiversity, communitymat)
Diversity.Ecology.generalisedpielou(metacommunityDiversity, communitymat)
```

!!! note
    `generalisedpielou` is not exported, so it must be qualified as above (or
    imported explicitly with
    `using Diversity.Ecology: generalisedpielou`). Every other measure on this
    page is exported by the submodule.

## Faith's PD

[`faith_pd`](@ref Diversity.Ecology.faith_pd) and
[`generalisedfaith_pd`](@ref Diversity.Ecology.generalisedfaith_pd) are the
exception on this
page, in two ways. The measures above are the classical ones for wholly distinct
types, and refuse to run when similarity is present; Faith's phylogenetic
diversity is the opposite - it is meaningless without a tree, so it is defined
only for a metacommunity built over the `PhyloBranches` that the `Phylo`
extension supplies, and therefore appears only once `Phylo` is loaded. See
[Phylogenetic diversity](phylogenetics.md).

It is also the one measure here that ignores abundances entirely: PD is the
total length of the branches spanned by the types that are *present*, so it
answers "how much evolutionary history is here?" rather than "how is it
distributed?".

```@contents
```

```@autodocs
Modules = [Diversity.Ecology]
Private = false
```

```@index
```
