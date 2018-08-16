In the **Diversity.Ecology** package, we replicate old ecological
diversity measures and generalised versions of them that relate to our
general measures of alpha, beta and gamma diversity at subcommunity
and ecosystem measures. The generalisations of the richness, Shannon
and Simpson are the only standard measures we are aware of whose
subcommunity components sum directly to the corresponding ecosystem
measure (although note that Simpson's index decreases for increased
diversity, so small components are more diverse).

#### Usage

Accessing the functionality in the package is simple:

```jldoctest
julia> using Diversity.Ecology, LinearAlgebra

julia> community = [10, 20, 20];

julia> community /= sum(community);

julia> diversity = simpson(community)
1×7 DataFrames.DataFrame
│ Row │ div_type │ measure │ type_level │ type_name │ partition_level │ partition_name │ diversity │
├─────┼──────────┼─────────┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Unique   │ Simpson │ types      │           │ subcommunity    │ 1              │ 0.36      │

julia> ecosystem = [2 2 0.; 0 2 2]';

julia> ecosystem /= sum(ecosystem);

julia> jaccard(ecosystem)
1×8 DataFrames.DataFrame
│ Row │ div_type │ measure │ q │ type_level │ type_name │ partition_level │ partition_name │ diversity │
├─────┼──────────┼─────────┼───┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Unique   │ Jaccard │ 0 │ types      │           │ metacommunity   │                │ 0.333333  │

julia> generalisedjaccard(ecosystem, [0, 1, 2])
3×8 DataFrames.DataFrame
│ Row │ div_type    │ measure │ q │ type_level │ type_name │ partition_level │ partition_name │ diversity │
├─────┼─────────────┼─────────┼───┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Arbitrary Z │ Jaccard │ 0 │ types      │           │ metacommunity   │                │ 0.333333  │
│ 2   │ Arbitrary Z │ Jaccard │ 1 │ types      │           │ metacommunity   │                │ 0.414214  │
│ 3   │ Arbitrary Z │ Jaccard │ 2 │ types      │           │ metacommunity   │                │ 0.5       │

julia> generalisedjaccard(ecosystem, [0, 1, 2], Matrix(1.0I, 3, 3))
3×8 DataFrames.DataFrame
│ Row │ div_type    │ measure │ q │ type_level │ type_name │ partition_level │ partition_name │ diversity │
├─────┼─────────────┼─────────┼───┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Arbitrary Z │ Jaccard │ 0 │ types      │           │ metacommunity   │                │ 0.333333  │
│ 2   │ Arbitrary Z │ Jaccard │ 1 │ types      │           │ metacommunity   │                │ 0.414214  │
│ 3   │ Arbitrary Z │ Jaccard │ 2 │ types      │           │ metacommunity   │                │ 0.5       │
```

```@contents
```

```@autodocs
Modules = [Diversity.Ecology]
Private = false
```

```@index
```
