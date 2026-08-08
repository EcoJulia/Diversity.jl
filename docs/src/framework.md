# The framework

This package is an executable form of a specific piece of mathematics. This page explains what that
mathematics measures, what each measure means, and where it came from — everything the rest of the
documentation assumes you already know.

## Origins

Three ideas, each generalising the one before it.

**Effective numbers** ([Hill, 1973](http://www.jstor.org/stable/1934352)). A diversity should be
reported as a *number of types*: a community of `S` equally abundant, completely distinct types has
diversity exactly `S`, whatever the order of the measure. That single requirement is what makes
diversities comparable, and it is the property everything here is built to preserve.

**Similarity** ([Leinster and Cobbold,
2012](http://www.esajournals.org/doi/abs/10.1890/10-2402.1)). Types are rarely wholly distinct. Two
species in the same genus, two age classes a year apart, two bacterial isolates differing by one
resistance — treating these as unrelated throws away most of what you know. Leinster and Cobbold
made diversity sensitive to a *similarity matrix* `Z`, while keeping the effective-number reading:
the answer becomes the effective number of *distinct* types.

**Partitioning** ([Reeve *et al.*](https://arxiv.org/abs/1404.6520)). A community is usually divided
— into sites, time points, host populations. Existing approaches aggregated across the whole
population, giving only "within", "between" and "total" diversity. This framework instead measures
each subcommunity *individually* and in the context of the whole, so subcommunities can be compared
directly with one another. That is what this package implements.

!!! note "Citing this work"
    Reeve R, Leinster T, Cobbold CA, Thompson J, Brummitt N, Mitchell SN, Matthews L.
    *How to partition diversity.* [arXiv:1404.6520](https://arxiv.org/abs/1404.6520).
    `CITATION.bib` in the repository has the BibTeX entry. The framework is independently implemented
    in R as [rdiversity](https://github.com/boydorr/rdiversity), against which this package is
    cross-validated.

## What you supply

A calculation needs three things, and they map directly onto the three arguments of
[`Metacommunity`](@ref):

- **types** — what the individuals are, and how similar they are to one another. That similarity is
  the matrix `Z`, where `Zᵢᵢ'` is the similarity between types `i` and `i'`, usually between 0 and 1.
  Types wholly distinct from one another means `Z = I`, which is what [`UniqueTypes`](@ref) provides;
- **a partition** — how the metacommunity divides into subcommunities ([`Subcommunities`](@ref), or
  [`Onecommunity`](@ref) for an undivided one);
- **abundances** — a matrix `P`, types down and subcommunities across.

⚠️ Abundances are **relative to the whole metacommunity** and sum to one across all of it, not one per
column. The column sums are then the subcommunity **weights** `w`, and the row sums the metacommunity
abundance `p`. Counts are normalised for you; floating-point abundances that miss are corrected with a
warning.

## Ordinariness

Everything is built from one derived quantity. The **ordinariness** of a type is `(Zp)ᵢ` — the
expected similarity between an individual of that type and an individual drawn at random from the
metacommunity. A type is ordinary when there is a lot of stuff like it about; its reciprocal is that
type's **uniqueness**.

```@repl framework
using Diversity
using Diversity.ShortNames
pop = [0.2 0.1; 0.1 0.3; 0.2 0.1]
getordinariness!(Metacommunity(pop))
```

With no similarity, ordinariness is just abundance — an individual is only like others of its own
type. Introduce similarity and each type inherits some of its neighbours' abundance:

```@repl framework
Z = [1.0 0.5 0.0; 0.5 1.0 0.5; 0.0 0.5 1.0]
getordinariness!(Metacommunity(pop, Z))
```

Every measure in the package is a ratio of ordinarinesses, averaged with a power mean. That is the
whole design.

## The viewpoint parameter

Every measure takes an order `q`, the **viewpoint parameter**, running from `0` to `∞`. It sets how
much weight rare types get:

| `q` | reads as | recovers |
|---|---|---|
| `0` | every type counts alike, however rare | richness |
| `1` | types counted in proportion to their abundance | Shannon entropy |
| `2` | commoner types weighted more heavily | Simpson's index |
| `∞` | only the commonest type matters | Berger–Parker |

```@repl framework
uneven = [0.5, 0.3, 0.15, 0.04, 0.01]
meta_gamma(Metacommunity(uneven), [0, 1, 2, Inf])[!, :diversity]
```

Five types are present, so at `q = 0` the diversity is 5; as `q` rises the two rare types stop
counting and it falls towards the effective number of common types. ⭐ `q = 1` is the usual default
when there is no reason to prefer another: it corresponds to Shannon entropy, the most studied case.

Diversities are **decreasing in `q`** — for α, ᾱ, ρ, ρ̄ and γ and their metacommunity counterparts.
⚠️ β and β̄ run the other way, *increasing* in `q`, being reciprocals of ρ and ρ̄; and the
metacommunity averages `B` and `B̄` are not monotone in either direction.

## The measures

Seven measures, each meaningful at two scales. ⭐ The two readings genuinely differ — a subcommunity's
gamma diversity is its *contribution* to the whole, while the metacommunity's gamma is the diversity
of the whole — so the table gives both:

| package | symbol | as a subcommunity measure | as a metacommunity measure |
|---|---|---|---|
| [`RawAlpha`](@ref) | α / A | estimate of naïve-community metacommunity diversity | naïve-community metacommunity diversity |
| [`NormalisedAlpha`](@ref) | ᾱ / Ā | diversity of the subcommunity in isolation | average diversity of the subcommunities |
| [`RawRho`](@ref) | ρ / R | **redundancy** of the subcommunity | average redundancy of the subcommunities |
| [`RawBeta`](@ref) | β / B | **distinctiveness** of the subcommunity | average distinctiveness of the subcommunities |
| [`NormalisedRho`](@ref) | ρ̄ / R̄ | **representativeness** of the subcommunity | average representativeness of the subcommunities |
| [`NormalisedBeta`](@ref) | β̄ / B̄ | estimate of effective number of distinct subcommunities | effective number of distinct subcommunities |
| [`Gamma`](@ref) | γ / G | contribution per individual toward metacommunity diversity | metacommunity similarity-sensitive diversity |

Two axes run through that table:

- **raw versus normalised** is a factor of the subcommunity's size `w`, and nothing else. A raw
  measure is *per individual*; a normalised one treats the subcommunity as a community in its own
  right.
- **ρ and β are reciprocals.** Redundancy asks how much of a subcommunity's diversity would survive
  its loss; distinctiveness asks how much is unique to it. They are the same question asked from
  opposite ends, which is why the package provides both.

`subdiv` and `metadiv` select the scale, and the wrapper functions name the pair directly —
`norm_sub_rho` is representativeness per subcommunity, `meta_gamma` is metacommunity diversity.

## What the measures show you

Each of the following is a published result of the framework, reproduced here by running it.

### Diversities are effective numbers

Five equally abundant, wholly distinct types have diversity exactly five — at every `q`. This is the
defining property, and everything else is built to preserve it.

```@repl framework
meta_gamma(Metacommunity(fill(1 / 5, 5)), [0, 1, 2, Inf])[!, :diversity]
```

### Representativeness is a proportion

Take a metacommunity where all types are equally abundant, and each subcommunity holds an equal share
of them. A subcommunity holding a fraction `r` of the types has representativeness exactly `r`:

```@repl framework
thirds = [1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1] ./ 6
norm_sub_rho(Metacommunity(thirds), 1)[!, :diversity]
```

Each of the three subcommunities holds two of the six types, and represents a third of the
metacommunity.

### Gamma is per individual, so it does not track alpha

⭐ This is the distinction that makes subcommunity gamma a new measure rather than a rescaling. Here
two subcommunities draw on types that are equally abundant in the metacommunity, but the first has
twice as many of them:

```@repl framework
sizes = [1 0; 1 0; 1 0; 1 0; 0 1; 0 1] ./ 6
norm_sub_alpha(Metacommunity(sizes), 1)[!, :diversity]
sub_gamma(Metacommunity(sizes), 1)[!, :diversity]
```

The first subcommunity is twice as diverse in isolation — but each of its individuals contributes
exactly as much to the metacommunity as each of the second's, because all six types are equally
abundant overall. Rarity, not richness, is what changes the contribution:

```@repl framework
rarity = [1 0 0; 1 0 0; 0 1 2; 0 1 2] ./ 8
norm_sub_alpha(Metacommunity(rarity), 1)[!, :diversity][1:2]
sub_gamma(Metacommunity(rarity), 1)[!, :diversity][1:2]
```

Now the first two subcommunities are the same size and equally diverse in isolation, but the first
one's types are three times rarer in the metacommunity — and its contribution is three times greater.

### Hidden diversity

⭐ The case the framework was built to expose. A subcommunity holding a single very rare type is
utterly dull in isolation, yet each of its individuals contributes enormously to the diversity of the
whole:

```@repl framework
rare = 1e-9
hidden = zeros(12, 3);
hidden[1, 1] = rare;               # one very rare type
hidden[2:11, 2] .= rare;           # ten equally rare types
hidden[12, 3] = 1 - 11rare;        # the common background
hid = Metacommunity(hidden);
norm_sub_alpha(hid, 1)[!, :diversity]
sub_gamma(hid, 1)[!, :diversity]
norm_sub_beta(hid, 1)[!, :diversity]
```

The first subcommunity has `ᾱ = 1` — one type, no diversity at all — but `γ ≈ 10⁹`. Adding nine more
equally rare types to the second subcommunity raises its alpha tenfold and leaves its gamma
untouched, because each individual is still just as rare; its distinctiveness falls, because a
ten-type subcommunity looks a little more like the metacommunity than a one-type one does. Alpha and
beta alone would have told you the first site was the least interesting in the study.

## Properties

**Invariance under shattering.** If you draw a boundary through a subcommunity that is internally
well mixed, you have not created a new subcommunity, and the answer should not change. The normalised
metacommunity measures — `Ā`, `R̄`, `B̄` and `G` — are invariant under such shattering:

```@repl framework
whole = [2 1; 1 3; 0 1] ./ 8
split = [2 0.5 0.5; 1 1.5 1.5; 0 0.5 0.5] ./ 8
norm_meta_beta(Metacommunity(whole), 1)[!, :diversity]
norm_meta_beta(Metacommunity(split), 1)[!, :diversity]
```

⚠️ The raw measures `A`, `R` and `B` are deliberately *not* invariant, and should not be: `R` is the
average redundancy of the subcommunities, and cutting one in two genuinely does create redundancy.

```@repl framework
raw_meta_rho(Metacommunity(whole), 1)[!, :diversity]
raw_meta_rho(Metacommunity(split), 1)[!, :diversity]
```

**Conditional independence.** A subcommunity's measures depend only on its own composition and on the
metacommunity as a whole — never on how the *rest* of the metacommunity happens to be divided up.
Without this, the apparent importance of one site would shift when someone redrew a boundary
elsewhere, and comparing subcommunities would be meaningless. It is what licenses the direct
comparisons above.

## Special cases

| when | what happens |
|---|---|
| `Z = I` (`UniqueTypes`) | the **naïve-type** case: types wholly distinct, and the measures reduce to Hill numbers |
| no shared types between subcommunities | the **naïve-community** case: `B = 1`, `B̄` is the effective number of subcommunities, and `G = A` |
| every subcommunity has the metacommunity's composition | **well-mixed**: `R̄ = B̄ = 1`, and `Ā = G` |
| one subcommunity | every beta measure is 1, and `α = ᾱ = γ` |

The naïve-type case recovers Hill numbers exactly:

```@repl framework
using Diversity.Hill
naive = [0.5, 0.3, 0.2]
meta_gamma(Metacommunity(naive), [0, 1, 2])[!, :diversity]
hillnumber(naive, [0, 1, 2])[!, :diversity]
```

### ⚠️ Similarity can break the bounds you expect

Read only the naïve-type case and you will absorb some inequalities that do not hold in general —
notably that representativeness cannot exceed 1. The paper's counterexample is a heterosexual
transmission network: six individuals as types, similarity being the chance of transmission, split
into a subcommunity of three males and one of three females. Each individual is *less* similar to its
own subcommunity than to the other one:

```@repl framework
Zsex = [1.0 0.5 0.5 0.7 0.7 0.7
        0.5 1.0 0.5 0.7 0.7 0.7
        0.5 0.5 1.0 0.7 0.7 0.7
        0.7 0.7 0.7 1.0 0.5 0.5
        0.7 0.7 0.7 0.5 1.0 0.5
        0.7 0.7 0.7 0.5 0.5 1.0]
sexes = [1 0; 1 0; 1 0; 0 1; 0 1; 0 1] ./ 6
sex = Metacommunity(sexes, Zsex);
norm_sub_rho(sex, 1)[!, :diversity]
norm_sub_beta(sex, 1)[!, :diversity]
```

Representativeness above 1, distinctiveness below it. ⚠️ Nor is this a pathological matrix — it
satisfies the triangle inequality — so no weaker assumption rescues those bounds. If your `Z` is not
the identity, check the behaviour you are relying on rather than assuming it.

```@contents
```

```@index
```
