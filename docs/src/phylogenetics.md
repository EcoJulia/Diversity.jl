In the **Phylogenetics** code, we generate Phylogenetic
diversity measures, based on Faith PD and extended by Chao.

#### Usage

Using the functionality in the package is simple:

 - Create a tree (using our Phylo package)
 - Create a PhyloTypes (AbstractTypes subtype) object from it
 - Create a Metacommunity from that
 - Calculate diversity!

```julia-repl
julia> using Phylo, Diversity

julia> species = ["Dog", "Cat", "Human", "Potato"];

julia> community = [4, 1, 3, 2] / 10;

julia> nt = rand(Nonultrametric(species))
NamedTree phylogenetic tree with 7 nodes and 6 branches
Leaf names:
String["Dog", "Cat", "Human", "Potato"]

julia> collect(getbranches(nt))
6-element Array{Pair{Int64,Phylo.Branch{String}},1}:
 [node "Node 2"]-->[0.43833500493556077 length branch 4]-->[node "Dog"]
 [node "Node 1"]-->[0.3488621822941236 length branch 2]-->[node "Cat"]
 [node "Node 2"]-->[0.46822991782974444 length branch 3]-->[node "Potato"]
 [node "Node 3"]-->[0.6558953328300103 length branch 5]-->[node "Node 1"]
 [node "Node 3"]-->[0.6667039144666622 length branch 6]-->[node "Node 2"]
 [node "Node 1"]-->[0.3165255563433085 length branch 1]-->[node "Human"]

julia> ph = PhyloTypes(nt);

julia> metaphylo = Metacommunity(community, ph);

julia> leafnames = gettypenames(metaphylo, true)
4-element Array{String,1}:
 "Dog"
 "Cat"
 "Human"
 "Potato"

julia> meta_gamma(metaphylo, 0)
1×7 DataFrames.DataFrame
│ Row │ measure │ q │ type_level │ type_name │ partition_level │
├─────┼─────────┼───┼────────────┼───────────┼─────────────────┤
│ 1   │ "Gamma" │ 0 │ "types"    │ ""        │ "metacommunity" │

│ Row │ partition_name │ diversity │
├─────┼────────────────┼───────────┤
│ 1   │ ""             │ 2.72761   │
```

```@contents
```

```@autodocs
Modules = [Diversity]
Private = false
```

```@index
```
