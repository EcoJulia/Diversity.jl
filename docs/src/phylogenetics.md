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
**processed** types — which is what the `raw` argument to the accessors selects:

```@repl phylo
gettypenames(metaphylo, true)
counttypes(metaphylo, true)
counttypes(metaphylo, false)
```

```@contents
```

```@index
```
