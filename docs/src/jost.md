Lou Jost's
[diversity](http://dx.doi.org/10.1111/j.2006.0030-1299.14714.x)
[measures](http://www.esajournals.org/doi/abs/10.1890/06-1736.1) are
found in the **Diversity.Jost** package.

#### Usage

Accessing the main functionality in the package is simple:

```jldoctest
julia> using Diversity.Jost

julia> ecosystem = [2 2 0; 0 2 2]';

julia> ecosystem /= sum(ecosystem);

julia> diversities = jostbeta(ecosystem, [0, 1, 2])
3×8 DataFrames.DataFrame
│ Row │ div_type │ measure  │ q │ type_level │ type_name │ partition_level │ partition_name │ diversity │
├─────┼──────────┼──────────┼───┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Unique   │ JostBeta │ 0 │ types      │           │ metacommunity   │                │ 1.5       │
│ 2   │ Unique   │ JostBeta │ 1 │ types      │           │ metacommunity   │                │ 1.41421   │
│ 3   │ Unique   │ JostBeta │ 2 │ types      │           │ metacommunity   │                │ 1.33333   │
```

```@contents
```

```@autodocs
Modules = [Diversity.Jost]
Private = false
```

```@index
```
