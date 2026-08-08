# Diversity.jl

## A package for measuring and partitioning diversity

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
functions exist both with unicode names (e.g. `ᾱ()`), which are
not automatically exported (as we feel they are too short) and with
matching longer ASCII names (e.g. `NormalisedAlpha()`), which are.
We also provide functions to calculate appropriate
`subdiv()` and `metadiv()` values for each measure, and a general
`diversity()` function to extract any diversity measure at a series of scales.

Accessing the main functionality in the package is simple:

```@repl usage
using Diversity
pop = [1 1 0; 2 0 0; 3 1 4]
pop = pop ./ sum(pop)
meta = Metacommunity(pop)
norm_meta_alpha(meta, [0, 1, 2, Inf])
```

Every measure returns a `DataFrame` in the same format, whatever the measure and
whatever the scale, so results can be compared and concatenated directly. The
`diversity` column holds the answer; the rest say what was calculated and for
what.

Adding a similarity matrix makes the measures similarity-sensitive — two types
that resemble each other now contribute less diversity between them than two
that do not:

```@repl usage
Z = [1.0 0 0; 0 1 1; 1 1 1]
meta_z = Metacommunity(pop, Z)
subdiv(RawRho(meta_z), 2)
```

Note that the abundances are relative to the **whole metacommunity** and must
sum to one across it — not one per subcommunity. Counts are normalised for you,
and floating point abundances that do not sum to one are corrected with a
warning.

```@contents
```

```@autodocs
Modules = [Diversity, Diversity.ShortNames]
Private = false
```

Private functions in module Diversity:

```@autodocs
Modules = [Diversity]
Public = false
```

```@index
```
