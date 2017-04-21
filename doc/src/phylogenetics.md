In the **Diversity.Phylogenetics** submodule, we generate Phylogenetic
diversity measures, based on Faith PD and extended by Chao.

#### Usage

Using the functionality in the package is simple:

 - Create a tree (using PhyloTrees)
 - Create a Phylogeny (AbstractTypes subtype) object from it
 - Create a Metacommunity from that

```
julia> using PhyloTrees

julia> using Diversity

julia> using Diversity.Phylogenetics

julia> species = ["Dog", "Cat", "Human"];

julia> community = [4, 3, 3];

julia> community /= sum(community);

julia> nt = NamedTree(species);

julia> n = addnode!(nt);

julia> addbranch!(nt, n, "Dog", 1.0);

julia> addbranch!(nt, n, "Cat", 1.0);

julia> r = addnode!(nt);

julia> addbranch!(nt, r, "Human", 2.0);

julia> addbranch!(nt, r, n, 1.0);

julia> nt
NamedTree phylogenetic tree with 5 nodes and 4 branches
Leaf names:
String["Human", "Cat", "Dog"]

julia> ph = Phylogeny(nt);

julia> leafnames = getphylonames(ph);

julia> order = mapreduce(name -> find(species .== name), append!, leafnames);

julia> metaphylo = Metacommunity(community[order], ph);

julia> meta_gamma(metaphylo, 0)
1×7 DataFrames.DataFrame
│ Row │ measure │ q │ type_level │ type_name │ partition_level │ partition_name │
├─────┼─────────┼───┼────────────┼───────────┼─────────────────┼────────────────┤
│ 1   │ "Gamma" │ 0 │ "types"    │ ""        │ "metacommunity" │ ""             │

│ Row │ diversity │
├─────┼───────────┤
│ 1   │ 2.5       │

```

```@contents
```

```@autodocs
Modules = [Diversity.Phylogenetics]
Private = false
```

```@index
```
