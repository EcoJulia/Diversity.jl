# Large metacommunities

This page is about measuring diversity over data that is large enough for the shape of the
calculation to matter — hundreds of thousands of subcommunities, thousands of types, or both. A
gridded landscape is the usual reason: a 1 km grid over a country is a few hundred thousand cells,
and a species list for it can run to thousands.

Nothing here is a different API. It is the same [`diversity`](@ref Diversity.diversity),
[`subdiv`](@ref Diversity.subdiv) and [`metadiv`](@ref Diversity.metadiv) as everywhere else, used
in the order that keeps the intermediate results small.

## What is actually big

Two things scale with the size of your data, and everything else is small beside them.

**Your abundance matrix**, which is types × subcommunities. That is yours, and it has to exist.

**The ordinariness**, `Z * p`, which is the same shape. When there is no similarity between types —
[`UniqueTypes`](@ref Diversity.UniqueTypes), the default — this is free, because the ordinariness
*is* the abundances and the package does not copy them. With a similarity matrix it is a second
array of the same size, computed once and cached on the metacommunity.

Everything else is proportional to the *answer* rather than to the data. The individual diversities
of a measure are computed on demand rather than stored, so building a measure over a metacommunity
you have already measured allocates nothing at all; a subcommunity-level result is one row per
subcommunity, at about 88 bytes a row.

So for `nt` types and `np` subcommunities, in double precision:

| | bytes |
|---|---|
| abundances | `8 nt np` |
| ordinariness, without similarity | 0 |
| ordinariness, with similarity | `8 nt np` |
| the similarity matrix itself | `8 nt²` |
| a `subdiv` result, per measure per order | `88 np` |
| an `inddiv` result, per measure per order | `88 nt np` |

A worked case: 2000 species on a 1 km grid over a box containing the UK is roughly 10⁶ cells, so
**16 GB** of abundances, or **32 GB** with a 2000 × 2000 similarity matrix — while a
subcommunity-level result for one measure at one order is about 88 MB.

## `inddiv` is the one that does not scale

[`inddiv`](@ref Diversity.inddiv) returns a row per type per subcommunity, not per subcommunity. At
the size above that is 2 × 10⁹ rows, which is not a table anyone wants. `subdiv` and `metadiv` are
the ones that scale; reach for `inddiv` on a subset, or on data small enough to look at.

## Build the metacommunity once

A [`Metacommunity`](@ref Diversity.Metacommunity) caches the three things every measure needs — the
ordinariness, the subcommunity weights and the metacommunity ordinariness — on first use. Each is a
full pass over the abundance matrix, so the first measure pays for them and every later one does
not.

The corollary is that a metacommunity is not a live view of your data. It computes from what it was
given and never invalidates, so build it when the abundances are final, and build a fresh one rather
than editing an existing one.

## Ask for everything at once

[`diversity`](@ref Diversity.diversity) takes a collection of levels, a collection of measures and
one or more orders, and returns them in a single table:

```@repl large
using Diversity, Diversity.ShortNames
pop = rand(20, 500); pop ./= sum(pop);
mc = Metacommunity(pop)
levels = [subcommunityDiversity, metacommunityDiversity];
result = diversity(levels, [ᾱ, ρ̄, Γ], mc, [0, 1, 2]);
size(result)
```

The arithmetic is the same as calling the wrappers one at a time — the power means still have to be
taken. What changes is how much is built on the way: one table rather than six that you then have to
join. Measured over 200 types and 200,000 subcommunities, for exactly that call:

| | time | memory |
|---|---|---|
| `diversity(levels, [ᾱ, ρ̄, Γ], mc, [0, 1, 2])` | 2035 ms | **165 MiB** |
| the six wrapper calls, then `vcat` | 2133 ms | 275 MiB |

## Do not build a table you are only going to write out

For a run whose output is destined for a file, a `DataFrame` is an expensive intermediate. Every
result function takes an optional first argument naming what to return, so the columns can be handed
to a writer without one:

```julia
using CSV, Tables
CSV.write("diversity.csv",
          subdiv(Tables.columntable, NormalisedAlpha(mc), [0, 1, 2]))
```

This is worth doing at scale because most of a result is repetition. Of the eight columns, five hold
one value repeated for every row and two cycle a short list of names; only `diversity` is data. They
are held as rules rather than arrays, and a column table keeps them that way, so `CSV.write` reads
each column once and never materialises it. Asking for a `DataFrame` — the default — materialises
them deliberately, so that what you get back is an ordinary mutable table.

!!! note "Which sinks get their own treatment"
    The first argument is resolved through `Tables.materializer`, and only types that define a
    method for it — `DataFrame` among them — are constructed directly. Anything else falls back to
    `Tables.columntable`, which is a perfectly good answer and is what you want here, but it does
    mean naming a writer rather than a table type will not do what it looks like it does. Write the
    file with `CSV.write` or `Arrow.write` and pass them a column table, as above.

## Empty subcommunities are nearly free

A subcommunity with no individuals in it — a sea cell in a species grid, a cell excluded from a
landscape simulation — has no defined diversity, and its rows come back as `NaN`:

```@repl large
sparse = zeros(4, 5);
sparse[:, [2, 4]] .= [1 4; 2 3; 3 2; 4 1];
sparse ./= sum(sparse);
subdiv(NormalisedAlpha(Metacommunity(sparse)), 1).diversity
```

Such subcommunities are recognised from the weights the metacommunity has already computed rather
than by inspecting their abundances, so they cost almost nothing: over 100,000 subcommunities of
which 27% were occupied, they accounted for 16% of the run before this was so, and now they account
for essentially none of it.

They also do not affect anything else. Removing them entirely gives *identical* answers for every
remaining subcommunity and for the metacommunity as a whole, at every order and with or without
similarity — abundances are relative to the whole metacommunity, so zeroes do not move the total; a
subcommunity's measures depend only on its own composition and on the metacommunity; and a power
mean ignores zero-weight entries. So you may keep a fixed grid and let cells fall in and out of use
without the measures noticing, or drop the empty cells and get the same numbers over a smaller
matrix — whichever suits the data you have.

```@contents
```

```@index
```
