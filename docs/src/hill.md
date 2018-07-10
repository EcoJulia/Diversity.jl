[Hill numbers](http://www.jstor.org/stable/1934352) are found in the
**Diversity.Hill** package.

#### Usage

Accessing the main functionality in the package is simple:

```jldoctest
ulia> using Diversity.Hill

julia> community = [10, 20, 20, 0, 3];

julia> community /= sum(community);

julia> diversities = hillnumber(community, [0, 1, 2])
3×8 DataFrames.DataFrame
│ Row │ div_type │ measure    │ q │ type_level │ type_name │ partition_level │ partition_name │ diversity │
├─────┼──────────┼────────────┼───┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Unique   │ HillNumber │ 0 │ types      │           │ subcommunity    │ 1              │ 4.0       │
│ 2   │ Unique   │ HillNumber │ 1 │ types      │           │ subcommunity    │ 1              │ 3.36264   │
│ 3   │ Unique   │ HillNumber │ 2 │ types      │           │ subcommunity    │ 1              │ 3.09021   │
```

```@contents
```

```@autodocs
Modules = [Diversity.Hill]
Private = false
```

```@index
```
