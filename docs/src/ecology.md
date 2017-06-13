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
julia> using Diversity.Ecology

julia> community = [10, 20, 20];

julia> community /= sum(community);

julia> diversity = simpson(community)
1×7 DataFrames.DataFrame
│ Row │ measure   │ q │ type_level │ type_name │ partition_level │
├─────┼───────────┼───┼────────────┼───────────┼─────────────────┤
│ 1   │ "Simpson" │ 2 │ "types"    │ ""        │ "subcommunity"  │

│ Row │ partition_name │ diversity │
├─────┼────────────────┼───────────┤
│ 1   │ "1"            │ 0.36      │

julia> ecosystem = [2 2 0.; 0 2 2]';

julia> ecosystem /= sum(ecosystem);

julia> jaccard(ecosystem)
1×7 DataFrames.DataFrame
│ Row │ measure   │ q │ type_level │ type_name │ partition_level │
├─────┼───────────┼───┼────────────┼───────────┼─────────────────┤
│ 1   │ "Jaccard" │ 0 │ "types"    │ ""        │ "metacommunity" │

│ Row │ partition_name │ diversity │
├─────┼────────────────┼───────────┤
│ 1   │ ""             │ 0.333333  │

julia> generalisedjaccard(ecosystem, [0, 1, 2])
3×7 DataFrames.DataFrame
│ Row │ measure   │ q │ type_level │ type_name │ partition_level │
├─────┼───────────┼───┼────────────┼───────────┼─────────────────┤
│ 1   │ "Jaccard" │ 0 │ "types"    │ ""        │ "metacommunity" │
│ 2   │ "Jaccard" │ 1 │ "types"    │ ""        │ "metacommunity" │
│ 3   │ "Jaccard" │ 2 │ "types"    │ ""        │ "metacommunity" │

│ Row │ partition_name │ diversity │
├─────┼────────────────┼───────────┤
│ 1   │ ""             │ 0.333333  │
│ 2   │ ""             │ 0.414214  │
│ 3   │ ""             │ 0.5       │

julia> generalisedjaccard(ecosystem, [0, 1, 2], eye(3))
3×7 DataFrames.DataFrame
│ Row │ measure   │ q │ type_level │ type_name │ partition_level │
├─────┼───────────┼───┼────────────┼───────────┼─────────────────┤
│ 1   │ "Jaccard" │ 0 │ "types"    │ ""        │ "metacommunity" │
│ 2   │ "Jaccard" │ 1 │ "types"    │ ""        │ "metacommunity" │
│ 3   │ "Jaccard" │ 2 │ "types"    │ ""        │ "metacommunity" │

│ Row │ partition_name │ diversity │
├─────┼────────────────┼───────────┤
│ 1   │ ""             │ 0.333333  │
│ 2   │ ""             │ 0.414214  │
│ 3   │ ""             │ 0.5       │
```

```@contents
```

```@autodocs
Modules = [Diversity.Ecology]
Private = false
```

```@index
```
