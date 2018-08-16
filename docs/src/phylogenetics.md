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
Creating Diversity to Phylo interface...

julia> species = ["Dog", "Cat", "Human", "Potato"];

julia> community = [4, 1, 3, 2] / 10;

julia> nt = rand(Nonultrametric(species))
BinaryTree{DataFrames.DataFrame,Dict{String,Any}} with 4 tips, 7 nodes and 6 branches.
Leaf names are Dog, Cat, Human and Potato


julia> collect(getbranches(nt))
6-element Array{Pair{Int64,Branch{String}},1}:
 [node "Node 2"]-->[0.5097599049872488 length branch 4]-->[node "Dog"]

 [node "Node 1"]-->[0.09388407179505037 length branch 2]-->[node "Human"]

 [node "Node 2"]-->[0.03694779049915409 length branch 3]-->[node "Potato"]

 [node "Node 3"]-->[1.1304368761257053 length branch 5]-->[node "Node 1"]

 [node "Node 3"]-->[2.690399241213393 length branch 6]-->[node "Node 2"]

 [node "Node 1"]-->[1.069077819992828 length branch 1]-->[node "Cat"]


julia> ph = PhyloTypes(nt);

julia> metaphylo = Metacommunity(community, ph);

julia> leafnames = gettypenames(metaphylo, true)
4-element Array{String,1}:
 "Dog"
 "Cat"
 "Human"
 "Potato"

julia> meta_gamma(metaphylo, 0)
1×8 DataFrames.DataFrame
│ Row │ div_type     │ measure │ q │ type_level │ type_name │ partition_level │ partition_name │ diversity │
├─────┼──────────────┼─────────┼───┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Phylogenetic │ Gamma   │ 0 │ types      │           │ metacommunity   │                │ 2.29217   │```

```@contents
```

```@autodocs
Modules = [Diversity]
Private = false
```

```@index
```
