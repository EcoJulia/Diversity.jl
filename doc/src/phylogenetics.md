In the **Diversity.Phylogenetics** submodule, we generate Phylogenetic
diversity measures, based on Faith PD and extended by Chao.

#### Usage

Using the functionality in the package is simple:

 - Create a tree (using our Phylo package)
 - Create a PhyloTypes (AbstractTypes subtype) object from it
 - Create a Metacommunity from that
 - Calculate diversity!

```
julia> using Phylo, Diversity, Diversity.Phylogenetics

julia> species = ["Dog", "Cat", "Human", "Potato"];

julia> community = [4, 1, 3, 2] / 10;

julia> nt = rand(Nonultrametric(species))
NamedTree phylogenetic tree with 5 nodes and 4 branches
Leaf names:
String["Dog", "Cat", "Human", "Potato"]

julia> collect(getbranches(nt))
4-element Array{Pair{Int64,Phylo.Branch{String}},1}:
 [node "Node 2"]-->[0.3546749125889703 length branch 4]-->[node "Cat"]
 [node "Node 1"]-->[0.46110670009439114 length branch 2]-->[node "Human"]
 [node "Node 2"]-->[0.20617420588508106 length branch 3]-->[node "Node 1"]
 [node "Node 1"]-->[0.7750408781162164 length branch 1]-->[node "Dog"]

julia> ph = PhyloTypes(nt);

julia> metaphylo = Metacommunity(community, ph);

julia> leafnames = gettypenames(metaphylo, true)
3-element Array{String,1}:
 "Dog"
 "Cat"
 "Human"
 
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
