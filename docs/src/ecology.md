# Diversity.Ecology

In the **Diversity.Ecology** submodule, we replicate old ecological
diversity measures and generalised versions of them that relate to our
general measures of alpha, beta and gamma diversity at subcommunity
and ecosystem measures. The generalisations of the richness, Shannon
and Simpson are the only standard measures we are aware of whose
subcommunity components sum directly to the corresponding ecosystem
measure (although note that Simpson's index decreases for increased
diversity, so small components are more diverse).

## Usage

Accessing the functionality in the package is simple:

```jldoctest
julia> using Diversity.Ecology, LinearAlgebra

julia> community = [10, 20, 20];

julia> community /= sum(community); #Convert to measurements to proportions

julia> diversity = simpson(community)
1×7 DataFrame
│ Row │ div_type │ measure │ type_level │ type_name │ partition_level │ partition_name │ diversity │
│     │ String   │ String  │ String     │ String    │ String          │ String         │ Float64   │
├─────┼──────────┼─────────┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Unique   │ Simpson │ types      │           │ subcommunity    │ 1              │ 0.36      │

julia> ecosystem = [2 2 0.; 0 2 2]';

julia> ecosystem /= sum(ecosystem);

julia> jaccard(ecosystem)
1×8 DataFrame
│ Row │ div_type │ measure │ q     │ type_level │ type_name │ partition_level │ partition_name │ diversity │
│     │ String   │ String  │ Int64 │ String     │ String    │ String          │ String         │ Float64   │
├─────┼──────────┼─────────┼───────┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Unique   │ Jaccard │ 0     │ types      │           │ metacommunity   │                │ 0.333333  │

julia> generalisedjaccard(ecosystem, [0, 1, 2])
3×8 DataFrame
│ Row │ div_type    │ measure │ q     │ type_level │ type_name │ partition_level │ partition_name │ diversity │
│     │ String      │ String  │ Int64 │ String     │ String    │ String          │ String         │ Float64   │
├─────┼─────────────┼─────────┼───────┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Arbitrary Z │ Jaccard │ 0     │ types      │           │ metacommunity   │                │ 0.333333  │
│ 2   │ Arbitrary Z │ Jaccard │ 1     │ types      │           │ metacommunity   │                │ 0.414214  │
│ 3   │ Arbitrary Z │ Jaccard │ 2     │ types      │           │ metacommunity   │                │ 0.5       │

julia> generalisedjaccard(ecosystem, [0, 1, 2], Matrix(1.0I, 3, 3))
3×8 DataFrame
│ Row │ div_type    │ measure │ q     │ type_level │ type_name │ partition_level │ partition_name │ diversity │
│     │ String      │ String  │ Int64 │ String     │ String    │ String          │ String         │ Float64   │
├─────┼─────────────┼─────────┼───────┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Arbitrary Z │ Jaccard │ 0     │ types      │           │ metacommunity   │                │ 0.333333  │
│ 2   │ Arbitrary Z │ Jaccard │ 1     │ types      │           │ metacommunity   │                │ 0.414214  │
│ 3   │ Arbitrary Z │ Jaccard │ 2     │ types      │           │ metacommunity   │                │ 0.5       │

julia> community = [0.7, 0.2, 0.1]

julia> pielou(community)
1×7 DataFrame
 Row │ div_type  measure  type_level  type_name  partition_level  partition_name  diversity 
     │ String    String   String      String     String           String          Float64   
─────┼──────────────────────────────────────────────────────────────────────────────────────
   1 │ Unique    Pielou   types                  subcommunity     1                0.729847

julia> communitymat = [10 20 30 20 0; #5 species (columns) and 6 sites (rows)
                       10 0 50 80 10;
                       60 10 90 0 0; 
                       10 10 10 10 10;
                       70 70 70 70 70;
                       10 0 0 90 0]'


```

```@contents
```

```@autodocs
Modules = [Diversity.Ecology]
Private = false
```

```@index
```
