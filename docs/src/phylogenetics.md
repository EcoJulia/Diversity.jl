# Phylogenetic diversity

When both `Diversity` and `Phylo` are loaded,
we generate Phylogenetic diversity measures, based on Faith PD
and extended by Chao.

## Usage

Using the functionality in the package is simple:

- Create a tree (using our Phylo package)
- Create a PhyloBranches (AbstractTypes subtype) object from it
- Create a Metacommunity from that
- Calculate diversity!

```@repl phylo
using Diversity, Phylo
species = ["Dog", "Human", "Cat"]
tree = RootedTree(species)
internal = createnode!(tree)
createbranch!(tree, internal, "Dog", 1.0)
createbranch!(tree, internal, "Human", 1.0)
root = createnode!(tree)
createbranch!(tree, root, internal, 1.0)
createbranch!(tree, root, "Cat", 2.0)
ph = PhyloBranches(tree)
metaphylo = Metacommunity([0.4, 0.3, 0.3], ph)
meta_gamma(metaphylo, 0)
```

The tree above is built explicitly so that the numbers on this page are stable;
`rand(Nonultrametric(species))` will give you a random tree to experiment with
instead.

## Raw and processed types

Diversity is measured over the *branches* of the tree rather than over the
species, since that is what carries the evolutionary history. The species you
supplied are still available as the **raw** types, and the branches are the
**processed** types - which is what the `raw` argument to the accessors selects:

```@repl phylo
gettypenames(metaphylo, true)
counttypes(metaphylo, true)
counttypes(metaphylo, false)
```

## Faith's PD

Faith's phylogenetic diversity is the older, historical measure: the **total
length of the branches** of the tree spanned by the species present, with no
normalisation. It lives in `Diversity.Ecology` with the other classical indices,
but works only on a metacommunity built over the [`PhyloBranches`](@ref) this
extension supplies - the calculation depends on how that type in particular maps
leaf abundances onto branches - and so appears only once `Phylo` is loaded:

```@repl phylo
using Diversity.Ecology: faith_pd, generalisedfaith_pd
generalisedfaith_pd(metacommunityDiversity, metaphylo)[!, :diversity]
```

The tree above has branches of length 1, 1, 1 and 2, so its total is 5.0.

There is no `q` to give it, and the output has no `q` column: Faith's PD is the
q = 0 case by definition, so there is no profile to ask for. It also does not
depend on the abundances at all, only on which species are present - which is
exactly what distinguishes it from the framework's own q = 0 diversity:

```@repl phylo
generalisedfaith_pd(metacommunityDiversity,
                    Metacommunity([0.98, 0.01, 0.01], ph))[!, :diversity]
meta_gamma(Metacommunity([0.98, 0.01, 0.01], ph), 0)[!, :diversity]
```

The first is unchanged, because the same three species are still there. The
second is not, because it measures diversity *per unit* of branch length - the
two differ by exactly the scale factor the phylogenetic types carry.

With several subcommunities, `faith_pd` gives the PD of each one in isolation:

```@repl phylo
faith_pd(Metacommunity([0.4 0.0; 0.1 0.2; 0.0 0.3], ph))[!,
                                                         [:partition_name,
                                                          :diversity]]
```

```@contents
```

```@index
```
