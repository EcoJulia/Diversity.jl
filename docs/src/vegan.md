# Coming from vegan

Many people arriving here already measure diversity with R's
[vegan](https://github.com/vegandevs/vegan). This page maps what you know onto what is here, and - as
importantly - says what is deliberately *absent* and why.

!!! note "How reliable is each row?"
    Rows marked ✅ are **checked by the test suite** on every run: `test/run_rcall.jl` computes both
    sides and asserts they agree. Rows marked ○ are our reading of vegan's documentation and have not
    been machine-checked - if one is wrong, please
    [tell us](https://github.com/EcoJulia/Diversity.jl/issues).

## Before the table: three differences that matter more

Read [What this package does differently](framework.md) first if you have not. In short:

1. **Beta diversity is not pairwise here.** There is no `dist` object to get back.
2. **Abundances are relative to the whole metacommunity**, not normalised per site.
3. **Ā * B̄ = G does not hold**, except at `q = 1`. There is no `adipart`/`multipart` equivalent, on
   purpose.

Everything below makes more sense once those are in place.

## Alpha diversity of single communities

| vegan | here | |
|---|---|---|
| `specnumber(x)` | `richness(x)` | ○ |
| `diversity(x, "shannon")` | `shannon(x)` - returns the *entropy*, as vegan does | ○ |
| `diversity(x, "simpson")` | `simpson(x)` | ○ |
| `diversity(x, "invsimpson")` | `meta_gamma(mc, 2)`, or `hillnumber(x, 2)` | ○ |
| `renyi(x, scales, hill = TRUE)` | `hillnumber(x, qs)`, or any measure over a vector of `q` | ○ |
| `renyi(x, scales)` | `log.(hillnumber(x, qs).diversity)` | ○ |
| `fisher.alpha`, `rarefy`, `specaccum` | No equivalent - this package measures diversity, it does not estimate unseen richness | |

```@repl vegan
using Diversity, Diversity.Ecology, Diversity.Hill
community = [10, 20, 20, 0, 3];
richness(community).diversity
shannon(community).diversity
hillnumber(community, [0, 1, 2]).diversity
```

`richness` counts the types actually present, so the zero above does not contribute - the answer is
4, not 5. That is `q = 0` behaving as it should.

## Dissimilarity between two communities

These are the genuinely pairwise measures, and they behave as vegan's do - but they take exactly two
subcommunities rather than returning a matrix over many.

| vegan | here | |
|---|---|---|
| `vegdist(x, "jaccard")` | `jaccard(x)` | ✅ |
| `vegdist(x, "gower")` (pre-2.7) | `gower(x, countzeros = true)` | ✅ |
| `vegdist(x, "altGower")` | `gower(x, countzeros = false)` | ✅ |
| `vegdist` with other methods | not provided | |

```@repl vegan
two = [2 2; 2 0; 0 2] ./ 8
jaccard(two).diversity
gower(two, countzeros = true).diversity
gower(two, countzeros = false).diversity
```

Note: **vegan 2.7 changed `method = "gower"`**, range-standardising columns first and dropping tied
columns from the denominator - which returns `NA` for two identical samples. This package keeps the
classic Gower (1971) reading, so `countzeros = true` matches *old* vegan. The cross-validation in
`test/run_rcall.jl` reproduces old vegan explicitly for this reason.

## Partitioning across many subcommunities

This is where the packages genuinely diverge, and where the extra capability is.

| vegan | here | |
|---|---|---|
| `adipart` (additive Ā + B̄ = G) | deliberately absent | |
| `multipart` (multiplicative Ā * B̄ = G) | deliberately absent except at `q = 1` | |
| `betadiver(x, method)` | no equivalent - pairwise beta | |
| `betadisper` | no equivalent | |
| - | `norm_sub_rho` - how *representative* each subcommunity is | |
| - | `raw_sub_beta` - how *distinctive* each subcommunity is | |
| - | `raw_sub_rho` - how *redundant* each subcommunity is | |
| - | `sub_gamma` - each subcommunity's *contribution* to the whole | |

The four rows with no vegan equivalent are the point of the package. Rather than one number for
"how much turnover is there overall", you get a value **per subcommunity**, comparable across
subcommunities, telling you which sites are distinctive, which are representative, and which
contribute most to the diversity of the whole. See
[Building a metacommunity](metacommunities.md) for a worked example that picks out each in turn.

```@repl vegan
sites = [10 0 0 5; 10 10 0 5; 0 10 10 5; 0 0 10 5]
mc = Metacommunity(sites)
norm_sub_rho(mc, 1).diversity     # the last site is the most representative
raw_sub_beta(mc, 1).diversity     # and the least distinctive
```

## Similarity between types

vegan treats species as wholly distinct, then handles functional or phylogenetic structure through
separate machinery. Here it is one argument.

| vegan | here | |
|---|---|---|
| `taxa2dist` + `taxondive` | a similarity matrix `Z`, via `GeneralTypes` | ○ |
| `treedive`, `treedist` | [Phylogenetic diversity](phylogenetics.md) with `PhyloBranches` | ○ |
| - | [Genetic diversity](genetics.md) from sequences or a VCF | |

Every measure in the package takes similarity, so there is no separate set of functions for
"functional diversity" or "phylogenetic diversity" - the same `norm_sub_rho` answers the taxonomic,
functional, phylogenetic and genetic version of the question depending only on the `Z` you supply.

## Data preparation

| vegan | here | |
|---|---|---|
| `decostand(x, "total")` | not needed - abundances are normalised on construction | |
| sites as **rows** | types as **rows**, subcommunities as **columns** - the transpose of vegan | |
| counts | passed directly; integers are normalised silently | |

Note: **The orientation is transposed relative to vegan**, which is the single most common early mistake.
vegan wants sites × species; a `Metacommunity` wants types × subcommunities. If your diversity values
look implausible, check that first.

```@contents
```

```@index
```
