# Diversity.jl

The main **Diversity** package provides basic numbers-equivalent
diversity measures (described in
[Hill, 1973](http://www.jstor.org/stable/1934352)),
similarity-sensitive diversity measures (generalised from Hill, and
described in
[Leinster and Cobbold, 2012](http://www.esajournals.org/doi/abs/10.1890/10-2402.1)),
and related alpha, beta and gamma diversity measures at the level of
the metacommunity and its component subcommunities (generalised in
turn from Leinster and Cobbold, and described in
[Reeve et al, 2014](http://arxiv.org/abs/1404.6520)). The diversity
functions exist both with unicode names (e.g. ```ᾱ()```), which are
not automatically exported (as we feel they are too short) and with
matching longer ASCII names (e.g. `NormalisedAlpha()`), which are.
We also provide functions to calculate appropriate
`subdiv()` and `metadiv()`
values for each measure, a general `diversity()` function for
extract any diversity measure at a series of scales.

Accessing the main functionality in the package is simple:

```jldoctest
julia> # Load the package into R
       using Diversity

julia> # Example population
       pop = [1 1 0; 2 0 0; 3 1 4];

julia> pop = pop / sum(pop);

julia> # Create Metacommunity object
       meta = Metacommunity(pop);

julia> diversities = norm_meta_alpha(meta, [0, 1, 2, Inf])
4×8 DataFrames.DataFrame
│ Row │ div_type │ measure         │ q   │ type_level │ type_name │ partition_level │ partition_name │ diversity │
├─────┼──────────┼─────────────────┼─────┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Unique   │ NormalisedAlpha │ 0.0 │ types      │           │ metacommunity   │                │ 2.16667   │
│ 2   │ Unique   │ NormalisedAlpha │ 1.0 │ types      │           │ metacommunity   │                │ 1.86121   │
│ 3   │ Unique   │ NormalisedAlpha │ 2.0 │ types      │           │ metacommunity   │                │ 1.63636   │
│ 4   │ Unique   │ NormalisedAlpha │ Inf │ types      │           │ metacommunity   │                │ 1.0       │

julia> Z = [1.0 0 0; 0 1 1; 1 1 1];

julia> meta_z = Metacommunity(pop, Z);

julia> rho = RawRho(meta_z);

julia> redundancies = subdiv(rho, 2)
3×8 DataFrames.DataFrame
│ Row │ div_type    │ measure │ q │ type_level │ type_name │ partition_level │ partition_name │ diversity │
├─────┼─────────────┼─────────┼───┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Arbitrary Z │ RawRho  │ 2 │ types      │           │ subcommunity    │ 1              │ 2.0       │
│ 2   │ Arbitrary Z │ RawRho  │ 2 │ types      │           │ subcommunity    │ 2              │ 3.0       │
│ 3   │ Arbitrary Z │ RawRho  │ 2 │ types      │           │ subcommunity    │ 3              │ 3.0       │
```

```@contents
```

```@autodocs
Modules = [Diversity]
Private = false
```

```@index
```
