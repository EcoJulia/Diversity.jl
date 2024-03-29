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

```julia-repl
julia> using Diversity.Ecology, LinearAlgebra

julia> community = [10, 20, 20];

julia> community = community ./ sum(community); # Convert counts to proportions

julia> diversity = simpson(community)
1×7 DataFrame
 Row │ div_type  measure  type_level  type_name  partition_level  partition_name  diversity 
     │ String    String   String      String     String           String          Float64   
─────┼──────────────────────────────────────────────────────────────────────────────────────
   1 │ Unique    Simpson  types                  subcommunity     1                    0.36

julia> ecosystem = [2 2 0; 0 2 2]';

julia> ecosystem = ecosystem ./ sum(ecosystem);

julia> jaccard(ecosystem)
1×8 DataFrame
 Row │ div_type  measure  q      type_level  type_name  partition_level  partition_name  diversity 
     │ String    String   Int64  String      String     String           String          Float64   
─────┼─────────────────────────────────────────────────────────────────────────────────────────────
   1 │ Unique    Jaccard      0  types                  metacommunity                     0.333333

julia> generalisedjaccard(Metacommunity(ecosystem))
1×7 DataFrame
 Row │ div_type  measure  type_level  type_name  partition_level  partition_name  diversity 
     │ String    String   String      String     String           String          Float64   
─────┼──────────────────────────────────────────────────────────────────────────────────────
   1 │ Unique    Jaccard  types                  metacommunity                     0.333333

julia> generalisedjaccard(ecosystem, Matrix(1.0I, 3, 3))
1×7 DataFrame
 Row │ div_type     measure  type_level  type_name  partition_level  partition_name  diversity 
     │ String       String   String      String     String           String          Float64   
─────┼─────────────────────────────────────────────────────────────────────────────────────────
   1 │ Arbitrary Z  Jaccard  types                  metacommunity                     0.333333

julia> community = [0.7, 0.2, 0.1];

julia> pielou(community)
1×7 DataFrame
 Row │ div_type  measure  type_level  type_name  partition_level  partition_name  diversity 
     │ String    String   String      String     String           String          Float64   
─────┼──────────────────────────────────────────────────────────────────────────────────────
   1 │ Unique    Pielou   types                  subcommunity     1                0.729847

julia> communitymat = [10 20 30 20 0; #5 sites/subcommunities (columns) and 6 species (rows)
                       10 0 50 80 10;
                       60 10 90 0 0; 
                       10 10 10 10 10;
                       70 70 70 70 70;
                       10 0 0 90 0];

julia> generalisedpielou(subcommunityDiversity, communitymat)
5×7 DataFrame
 Row │ div_type     measure  type_level  type_name  partition_level  partition_name  diversity 
     │ String       String   String      String     String           String          Float64   
─────┼─────────────────────────────────────────────────────────────────────────────────────────
   1 │ Arbitrary Z  Pielou   types                  subcommunity     1                0.781115
   2 │ Arbitrary Z  Pielou   types                  subcommunity     2                0.745557
   3 │ Arbitrary Z  Pielou   types                  subcommunity     3                0.888073
   4 │ Arbitrary Z  Pielou   types                  subcommunity     4                0.864562
   5 │ Arbitrary Z  Pielou   types                  subcommunity     5                0.622366

julia> generalisedpielou(metacommunityDiversity, communitymat)
1×7 DataFrame
 Row │ div_type     measure  type_level  type_name  partition_level  partition_name  diversity 
     │ String       String   String      String     String           String          Float64   
─────┼─────────────────────────────────────────────────────────────────────────────────────────
   1 │ Arbitrary Z  Pielou   types                  metacommunity                     0.510146
```

```@contents
```

```@autodocs
Modules = [Diversity.Ecology]
Private = false
```

```@index
```
