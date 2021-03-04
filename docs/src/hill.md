# Diversity.Hill

[Hill numbers](http://www.jstor.org/stable/1934352) are found in the
**Diversity.Hill** submodule.

#### Usage

Accessing the main functionality in the package is simple:

```jldoctest
julia> using Diversity.Hill

julia> community = [10, 20, 20, 0, 3];

julia> community /= sum(community);

julia> diversities = hillnumber(community, [0, 1, 2])
3×7 DataFrame
│ Row │ measure    │ q     │ type_level │ type_name │ partition_level │ partition_name │ diversity │
│     │ String     │ Int64 │ String     │ String    │ String          │ String         │ Float64   │
├─────┼────────────┼───────┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ HillNumber │ 0     │ types      │           │ subcommunity    │ 1              │ 4.0       │
│ 2   │ HillNumber │ 1     │ types      │           │ subcommunity    │ 1              │ 3.36264   │
│ 3   │ HillNumber │ 2     │ types      │           │ subcommunity    │ 1              │ 3.09021   │
```

```@contents
```

```@autodocs
Modules = [Diversity.Hill]
Private = false
```

```@index
```
