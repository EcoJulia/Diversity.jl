# Building a metacommunity

Everything this package computes starts from a `Metacommunity`, and it is made of three things:

- **types** — what the individuals are, and how similar they are to each other;
- **a partition** — how the whole is divided into subcommunities;
- **abundances** — how much of each type is in each subcommunity.

If you supply only abundances, the other two are inferred: types become
[`UniqueTypes`](@ref) (all wholly distinct from one another), and the partition becomes one
subcommunity per column.

## You have X, write Y

```@repl building
using Diversity
using LinearAlgebra
counts = [10 0 0; 5 5 0; 0 5 10]          # 3 types down, 3 subcommunities across
Z = [1.0 0.5 0.0; 0.5 1.0 0.5; 0.0 0.5 1.0]
nothing # hide
```

| what you have | what to write | what you get |
|---|---|---|
| one community, no similarity | `Metacommunity(vec)` | `UniqueTypes`, `Onecommunity` |
| several subcommunities, no similarity | `Metacommunity(matrix)` | `UniqueTypes`, one subcommunity per column |
| a similarity matrix as well | `Metacommunity(matrix, Z)` | `GeneralTypes` built from `Z` |
| named types or subcommunities | `Metacommunity(matrix, types, partition)` | exactly what you passed |
| a phylogeny, sequences or a VCF | `Metacommunity(matrix, PhyloBranches(tree))` etc. | see [Phylogenetic](phylogenetics.md) and [Genetic diversity](genetics.md) |
| an `EcoBase` assemblage | `Metacommunity(assemblage)` | its occurrences and places, with `UniqueTypes` — or `GeneralTypes` holding its similarity, if it has any |
| the same structure, new numbers | `Metacommunity(newabundances, oldmeta)` | reuses the types and partition |

```@repl building
Metacommunity(counts)
Metacommunity(counts, Z)
Metacommunity(counts, UniqueTypes(["ash", "oak", "elm"]),
              Subcommunities(["north", "middle", "south"]))
```

Naming your types and subcommunities costs one line and pays for itself the moment you read a
results table, because the names appear in the `type_name` and `partition_name` columns.

## Counts or proportions

You can pass either, but they are treated differently, and the difference is deliberate.

```@repl building
getabundance(Metacommunity(counts))          # counts: normalised silently
```

Integer counts are counts, so they are divided by their total without comment. Floating-point
abundances are assumed to be relative already, so if they do not sum to one you get a warning:

```@repl building
percolumn = counts ./ sum(counts, dims = 1)   # each column sums to 1 - a common mistake
sum(percolumn)                                # ... so the whole thing sums to 3, not 1
```

Handing that to `Metacommunity` gets you a warning and a correction:

```julia
julia> Metacommunity(percolumn);
┌ Warning: Abundances not normalised to 1, correcting...
└ @ Diversity Metacommunity.jl:124
```

Abundances here are relative to the
**whole metacommunity**, not to each subcommunity — the column sums are the subcommunity *weights*,
and normalising per column throws them away, making every subcommunity look the same size. If you
see it, you should probably go back to the original data and normalise over the whole metacommunity.

## A worked example

Ten sites and nine species, built so you can predict what each measure will say:

- **Sites 1–8** lie along a gradient, each holding two or three species, turning over as you move
  along it.
- **Reference** holds all eight gradient species, evenly.
- **Refuge** is tiny and holds a single species found nowhere else.

```@repl building
species = ["Sp $c" for c in "ABCDEFGH"];
push!(species, "Refuge sp");
sites = ["Site $i" for i in 1:8];
push!(sites, "Reference"); push!(sites, "Refuge");
abundance = zeros(Int, 9, 10);
for j in 1:8, i in max(1, j - 1):min(8, j + 1)
    abundance[i, j] = 10                      # a moving window along the gradient
end
abundance[1:8, 9] .= 6;                       # Reference: every species, evenly
abundance[9, 10] = 4;                         # Refuge: one species, found nowhere else
mc = Metacommunity(abundance, UniqueTypes(species), Subcommunities(sites))
```

That structure is much easier to see than to read off the loop:

```@example building
using Plots
heatmap(sites, species, abundance, yflip = true, xrotation = 45, c = :Blues,
        title = "Individuals per species per site", colorbar_title = "count")
```

The dark band down the diagonal is the gradient: each site shares species with its neighbours and
none at all with the far end. *Reference* is the one column touching every gradient species, at lower
abundance. *Refuge* is the single cell in a row of its own — nothing else is in that site, and that
species is nowhere else.

Before running anything, decide what you expect. *Reference* should be the most diverse in isolation
and the most representative; *Refuge* should be the least diverse in isolation but the most
distinctive, and should contribute most per individual. Now check:

```@repl building
using DataFrames
results = DataFrame(site = sites,
                    weight = round.(getweight(mc), digits = 3),
                    ᾱ = round.(norm_sub_alpha(mc, 1).diversity, digits = 2),
                    ρ̄ = round.(norm_sub_rho(mc, 1).diversity, digits = 2),
                    β = round.(raw_sub_beta(mc, 1).diversity, digits = 2),
                    γ = round.(sub_gamma(mc, 1).diversity, digits = 1))
```

Two things worth reading off that table:

- **Refuge scores `β = 1`, the maximum possible.** Distinctiveness is 1 exactly when nothing outside
  the subcommunity resembles anything inside it — which is true here by construction. Its
  representativeness is correspondingly at *its* minimum, which is the subcommunity's own weight.
- **Refuge has the lowest `ᾱ` and much the highest `γ`.** A single-species site is as dull as a
  site can be in isolation, yet each of its individuals contributes far more to the diversity of the
  whole than an individual from anywhere else. Alpha and beta alone would have told you to ignore it.

## How the partition changes the answer

The same abundances, divided three ways.

**Undivided.** With one subcommunity there is nothing to be distinct *from*, so every beta measure
collapses to 1 and alpha equals gamma:

```@repl building
undivided = Metacommunity(vec(sum(abundance, dims = 2)), UniqueTypes(species),
                          Onecommunity())
norm_sub_beta(undivided, 1).diversity
norm_sub_rho(undivided, 1).diversity
```

**Divided.** The gamma diversity of the whole is unchanged — dividing a community does not alter what
is in it — but the beta measures now carry information:

```@repl building
meta_gamma(undivided, 1).diversity
meta_gamma(mc, 1).diversity
norm_meta_beta(mc, 1).diversity
```

**Shattered.** Now split *Reference* into two identical halves. This creates no new ecology: the two
halves have the same composition as their parent, so an honest measure of "how many distinct
subcommunities are there?" must not move.

```@repl building
halves = hcat(abundance[:, 1:8], abundance[:, 9] .÷ 2, abundance[:, 9] .÷ 2,
              abundance[:, 10])
shattered = Metacommunity(halves, UniqueTypes(species), Subcommunities(11))
norm_meta_beta(mc, 1).diversity, norm_meta_beta(shattered, 1).diversity
norm_meta_alpha(mc, 1).diversity, norm_meta_alpha(shattered, 1).diversity
```

The normalised measures do not move. The raw ones deliberately do:

```@repl building
raw_meta_rho(mc, 1).diversity, raw_meta_rho(shattered, 1).diversity
```

**This is the clearest way to see what "raw" and "normalised" mean.** Raw redundancy *should* rise
when you cut a subcommunity in two, because you have genuinely created two subcommunities that
duplicate each other. Normalised measures answer the question "how many *distinct* subcommunities are
there really?", and the answer is unchanged. The two differ by exactly the subcommunity's weight `w`:
a raw measure keeps it and so sees size, a normalised one divides it out and cannot. Which you want depends on the question you are asking;
see [What this package does differently](framework.md) for why the framework guarantees the second.

## Getting the results out

Every measure returns a long `DataFrame` in the same shape, whatever the measure and whatever the
scale. That is convenient for combining results and inconvenient for reading, so to get back to a
familiar types × subcommunities matrix, `unstack` it:

```@repl building
ind = inddiv(NormalisedAlpha(mc), 1)
unstack(ind, :type_name, :partition_name, :diversity)
```

The columns that identify a row are `div_type`, `measure`, `q`, `type_level`, `type_name`,
`partition_level` and `partition_name`; `diversity` holds the answer. Selecting several orders at
once and filtering afterwards is usually faster than repeated calls, because the individual
diversities are computed once when the measure is built:

```@repl building
profile = subdiv(NormalisedAlpha(mc), [0, 1, 2]);
filter(row -> row.partition_name == "Reference", profile)
```

## Now on real data

The dataset above was built so you could check the measures against your own expectations. Real data
is the other way round — you do not know the answer, which is the point of measuring. Here is the
same analysis on the European amphibian distributions that ship with
[SpatialEcology](https://github.com/EcoJulia/SpatialEcology.jl): 73 species across 1010 grid cells.

Any `EcoBase` assemblage can be measured directly, with no conversion:

```@repl amphibian
using Diversity, SpatialEcology, CSV, DataFrames
file = joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv");
raw = CSV.read(file, DataFrame);
amph = Assemblage(raw[!, 4:end], raw[!, 1:3], sitecolumns = false)
countsubcommunities(amph), counttypes(amph)
```

```@repl amphibian
rho = norm_sub_rho(amph, 1)[!, :diversity];
gamma = sub_gamma(amph, 1)[!, :diversity];
extrema(rho)
extrema(gamma)
```

Representativeness runs from 0.0015 to 0.76: some cells hold an assemblage much like Europe's as a
whole, others almost nothing like it. Contribution per individual spans two orders of magnitude —
the cells at the top are those holding species found almost nowhere else, which is exactly the
"hidden diversity" the Refuge site showed in miniature.

Because the assemblage carries its own coordinates, the results can go straight onto a map:

```@example amphibian
using Plots
plot(norm_sub_rho(amph, 1), amph, title = "Representativeness of European amphibian faunas",
     markersize = 2)
```

Read that as a question about *reserve selection*: the dark cells are the ones whose amphibian
fauna is least like Europe's overall. They are not necessarily the most species-rich — that is what
`norm_sub_alpha` would show — and the difference between the two maps is precisely what the framework
was built to expose.

```@contents
```

```@index
```
