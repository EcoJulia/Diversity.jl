# Diversity.Phylogenetics

In the **Diversity** module's **Diversity.Phylogenetics** code,
we generate Phylogenetic diversity measures, based on Faith PD
and extended by Chao.

## Usage

Using the functionality in the package is simple:

- Create a tree (using our Phylo package)
- Create a PhyloTypes (AbstractTypes subtype) object from it
- Create a Metacommunity from that
- Calculate diversity!

```julia-repl
julia> using Diversity, Phylo, Diversity.Phylogenetics

julia> species = ["Dog", "Cat", "Human", "Potato"];

julia> community = [4, 1, 3, 2] / 10;

julia> nt = rand(Nonultrametric(species))
RootedTree with 4 tips, 7 nodes and 6 branches.
Leaf names are Dog, Human, Potato and Cat

julia> collect(getbranches(nt))
6-element Vector{LinkBranch{OneRoot, String, Dict{String, Any}, Float64}}:
 LinkBranch 7, from node Node 5 to node Human (length 0.48129057007442144).

 LinkBranch 8, from node Node 5 to node Potato (length 0.04547304399403693).

 LinkBranch 9, from node Node 6 to node Node 5 (length 0.5257296693452039).

 LinkBranch 10, from node Node 6 to node Cat (length 1.2058770065006985).

 LinkBranch 11, from node Node 7 to node Dog (length 1.0055430448967724).

 LinkBranch 12, from node Node 7 to node Node 6 (length 0.07011899861455448).


julia> ph = PhyloBranches(nt);

julia> metaphylo = Metacommunity(community, ph);

julia> leafnames = gettypenames(metaphylo, true)
4-element Vector{String}:
 "Dog"
 "Human"
 "Potato"
 "Cat"

julia> meta_gamma(metaphylo, 0)
1×8 DataFrame
│ Row │ div_type            │ measure │ q     │ type_level │ type_name │ partition_level │ partition_name │ diversity │
│     │ String              │ String  │ Int64 │ String     │ String    │ String          │ String         │ Float64   │
├─────┼─────────────────────┼─────────┼───────┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Phylogenetic Branch │ Gamma   │ 0     │ types      │           │ metacommunity   │                │ 3.48192   │
```

```@contents
```

```@autodocs
Modules = [Diversity.Phylogenetics]
Private = false
```

```@index
```
