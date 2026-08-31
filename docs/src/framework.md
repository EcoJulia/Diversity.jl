# The framework

This package is an executable form of a specific piece of mathematics. This page explains what that
mathematics measures, what each measure means, and where it came from - everything the rest of the
documentation assumes you already know.

## Origins

Three ideas, each generalising the one before it.

**Effective numbers** ([Hill, 1973](http://www.jstor.org/stable/1934352)). A diversity should be
reported as a *number of types*: a community of `S` equally abundant, completely distinct types has
diversity exactly `S`, whatever the order of the measure. That single requirement is what makes
diversities comparable, and it is the property everything here is built to preserve. Hill did not
invent the underlying mathematics - he took a family of entropy measures from information theory and
reported them in units an ecologist can use, rather than in bits. That is why exponentiating turns up
so often here: a Hill number is the exponential of an entropy. Where that family came from is at the
end of this section.

**Similarity** ([Leinster and Cobbold,
2012](http://www.esajournals.org/doi/abs/10.1890/10-2402.1)). Types are rarely wholly distinct. Two
species in the same genus, two age classes a year apart, two bacterial isolates differing by one
resistance - treating these as unrelated throws away most of what you know. Leinster and Cobbold
generalised Hill's mathematics to make diversity sensitive to a *similarity matrix* `Z`, while
keeping the effective-number reading: the answer becomes the effective number of *distinct* types.

**Partitioning** ([Reeve *et al.*](https://arxiv.org/abs/1404.6520)). A community is usually
divided - into sites, time points, host populations. Existing approaches aggregated across the whole
population, giving only "within", "between" and "total" diversity. This framework generalises
Leinster and Cobbold's mathematics in turn, measuring each subcommunity *individually* and in the
context of the whole, so subcommunities can be compared directly with one another. That is what this
package implements.

**Where the mathematics came from.**
[Shannon (1948)](https://doi.org/10.1002/j.1538-7305.1948.tb01338.x) defined the entropy of a
probability distribution, to measure information in a communication channel.
[Rényi (1961)](https://projecteuclid.org/euclid.bsmsp/1200512181) generalised it, showing that
Shannon's entropy is one member of a whole *family* indexed by a parameter - the `q` that runs
through everything above.

!!! note "Citing this work"
    Reeve R, Leinster T, Cobbold CA, Thompson J, Brummitt N, Mitchell SN, Matthews L.
    *How to partition diversity.* [arXiv:1404.6520](https://arxiv.org/abs/1404.6520).
    `CITATION.bib` in the repository has the BibTeX entry. The framework is independently implemented
    in R as [rdiversity](https://github.com/boydorr/rdiversity), against which this package is
    cross-validated.

## What this package does differently

If you are arriving from a classical diversity package - R's
[vegan](https://github.com/vegandevs/vegan) especially - three things here will not be what you
expect. None of them is an accident, and knowing them up front removes most of the surprise.

**Beta diversity is not pairwise.** `vegdist` and `betadiver` give you a dissimilarity *between two
sites*, and a matrix of them for a whole dataset. Here, every beta measure compares one subcommunity
*with the metacommunity it belongs to*, giving one value per subcommunity rather than a matrix. That
is what makes subcommunities directly comparable with one another - and it is why
[`jaccard`](ecology.md), which genuinely is a pairwise index, only accepts exactly two subcommunities.

**Abundances are relative to the whole metacommunity.** A classical package takes a sites × species
table and normalises each row. Here the whole matrix sums to one, so the column sums are the
subcommunity **weights**, and a subcommunity's size is part of the answer rather than something
divided out. Give it counts and it will normalise them for you; give it proportions that sum to one
*per column* and you will get a warning telling you they did not sum to one overall.

Note: **α * β = γ does not hold here, and that is the central design decision.** vegan offers
`adipart` and `multipart` for exactly that partitioning, and a great deal of the literature assumes
it. In this framework the relationship holds when `q = 1`, and in degenerate cases, but not in
general.

What was bought by giving it up:

| property | what it means |
|---|---|
| invariance under shattering | splitting a well-mixed subcommunity in two does not change the normalised metacommunity measures - an arbitrary boundary cannot manufacture diversity |
| conditional independence | a subcommunity's measures depend only on itself and on the metacommunity as a whole, never on how the *rest* is divided up |
| comparability | consequently, two subcommunities can be compared directly, and ranked |

Those three are what let you ask *"which of my sites is the most distinctive?"* and get an answer
that does not change when someone redraws a boundary elsewhere. The conflict is with the first of
them: solving the partition for beta makes it gamma divided by alpha, which ties its behaviour under
shattering to alpha's. Gamma cannot depend on how the metacommunity was divided, so an alpha that is
not invariant forces a beta that is not invariant either, and the size-weighted average alpha the
partition is usually built on is exactly that. Reeve *et al.* prove this for Jost's measures, by
splitting one well-mixed subcommunity of a two-subcommunity metacommunity: their alpha and beta both
move for every `q` except 1, where they coincide with the measures used here. See
[Properties](@ref) below for what this framework guarantees instead.

## What you supply

A calculation needs three things, and they map directly onto the three arguments of
[`Metacommunity`](@ref):

- **types** - what the individuals are, and how similar they are to one another. That similarity is
  the matrix `Z`, where `Zᵢⱼ` is the similarity between types `i` and `j`, usually between 0 and 1.
  Types wholly distinct from one another means `Z = I`, which is what [`UniqueTypes`](@ref) provides;
- **a partition** - how the metacommunity divides into subcommunities ([`Subcommunities`](@ref), or
  [`Onecommunity`](@ref) for an undivided one);
- **abundances** - a matrix `P`, types down and subcommunities across.

Note:Abundances are **relative to the whole metacommunity** and sum to one across all of it, not one per
column. The column sums are then the subcommunity **weights** `w`, and the row sums the metacommunity
abundance `p`. Counts are normalised for you; floating-point abundances that miss are corrected with a
warning.

## Ordinariness

Everything is built from one derived quantity. The **ordinariness** of a type is `(Zp)ᵢ` - the
expected similarity between an individual of that type and an individual drawn at random from the
metacommunity. A type is ordinary when there is a lot of stuff like it about; its reciprocal is that
type's **uniqueness**.

```@repl framework
using Diversity
using Diversity.ShortNames
pop = [0.2 0.1; 0.1 0.3; 0.2 0.1]
getordinariness!(Metacommunity(pop))
```

With no similarity, ordinariness is just abundance - an individual is only like others of its own
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
| `∞` | only the commonest type matters | Berger-Parker |

```@repl framework
uneven = [0.5, 0.3, 0.15, 0.04, 0.01]
meta_gamma(Metacommunity(uneven), [0, 1, 2, Inf])[!, :diversity]
```

Five types are present, so at `q = 0` the diversity is 5; as `q` rises the two rare types stop
counting and it falls towards the effective number of common types. `q = 1` is the usual default
when there is no reason to prefer another: it corresponds to Shannon entropy, the most studied case.

Diversities are **decreasing in `q`** - for α, ᾱ, ρ, ρ̄ and γ and their metacommunity counterparts.
Note: β and β̄ run the other way, *increasing* in `q`, being reciprocals of ρ and ρ̄; and the
metacommunity averages `B` and `B̄` are not monotone in either direction.

### What "viewpoint" means

The name is the useful part. `q` is not a property of the community; it is a property of *the person
asking*. A conservationist counting how many species are present at all takes `q = 0`, giving a rare
species and a dominant one equal say. Someone asking what an ecosystem is functionally *made of*
takes a high `q`, where a species at 0.1% abundance barely registers. Neither is more correct - they
are different viewpoints on the same data, and the framework makes you say which you are taking.

That is why a single number is usually the wrong output. Report a **diversity profile** instead:

```@example profile
using Diversity, Plots
uneven = [0.5, 0.3, 0.15, 0.04, 0.01]
even = fill(0.2, 5)
qs = 0:0.1:5
plot(qs, meta_gamma(Metacommunity(uneven), qs)[!, :diversity],
     label = "uneven", xlabel = "viewpoint parameter q",
     ylabel = "diversity", linewidth = 2)
plot!(qs, meta_gamma(Metacommunity(even), qs)[!, :diversity],
      label = "even", linewidth = 2)
```

Both communities contain five types, so both start at 5. The even one stays there - with nothing rare
to discount, every viewpoint agrees. The uneven one falls away as `q` rises, and *how fast it falls
is the evenness*. A profile therefore says more than richness and evenness reported separately,
because it shows them as one curve.

**If you know vegan, you have met this already**: `renyi(x, scales, hill = TRUE)` computes exactly
this, and its `scales` argument is `q`. `renyiaccum` plots the profile. The difference here is that
you can take a profile of *any* measure, not just gamma - a profile of `norm_sub_rho` shows how a
subcommunity's representativeness depends on whether you care about its rare types.

**When in doubt use `q = 1`.** It is the only order weighting each type exactly in proportion to
its abundance, it corresponds to Shannon entropy, and it is where the multiplicative relationships
hold. Every example in this documentation uses it unless there is a reason not to.

## The measures

Seven measures, each meaningful at two scales. The two readings genuinely differ - a subcommunity's
gamma diversity is its *contribution* to the whole, while the metacommunity's gamma is the diversity
of the whole - so the table gives both:

| package | symbol | as a subcommunity measure | as a metacommunity measure |
|---|---|---|---|
| [`RawAlpha`](@ref) | α / A | estimate of naive-community metacommunity diversity | naive-community metacommunity diversity |
| [`NormalisedAlpha`](@ref) | ᾱ / Ā | diversity of the subcommunity in isolation | average diversity of the subcommunities |
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

`subdiv` and `metadiv` select the scale, and the wrapper functions name the pair directly -
`norm_sub_rho` is representativeness per subcommunity, `meta_gamma` is metacommunity diversity.

## What the measures show you

Each of the following is a published result of the framework, reproduced here by running it.

### Diversities are effective numbers

Five equally abundant, wholly distinct types have diversity exactly five - at every `q`. This is the
defining property, and everything else is built to preserve it.

```@repl framework
meta_gamma(Metacommunity(fill(1 / 5, 5)), [0, 1, 2, Inf])[!, :diversity]
```

### Representativeness is a proportion

Take a metacommunity where all types are equally abundant, and each subcommunity holds an equal share
of them. A subcommunity holding a fraction `r` of the types has representativeness exactly `r`:

```@repl framework
twothirds = [1 1 0; 1 1 0; 0 1 1; 0 1 1; 1 0 1; 1 0 1] ./ 12
norm_sub_rho(Metacommunity(twothirds), 1)[!, :diversity]
```

Each type sits in two of the three subcommunities, so each subcommunity holds four of the six types -
two thirds of them - and represents two thirds of the metacommunity. Note: It is still only a *third* of
the individuals, and the two numbers are deliberately different: it is the share of the types that
the answer follows, not the share of the population.

Size is irrelevant in the other direction too. Give every subcommunity every type, in the
metacommunity's own proportions, and representativeness is 1 however uneven their sizes are:

```@repl framework
allshared = [1 2 3; 1 2 3; 1 2 3; 1 2 3; 1 2 3; 1 2 3] ./ 36
getweight(Metacommunity(allshared))
norm_sub_rho(Metacommunity(allshared), 1)[!, :diversity]
```

The last subcommunity has three times the individuals of the first and both are perfectly
representative. A subcommunity is representative when it *looks like* the metacommunity, not when
it is a large part of it.

**That insensitivity is what "normalised" means.** Representativeness is
the normalised measure of the pair: the subcommunity's weight `w` is divided out, so it is treated as
a community in its own right and its size cannot reach the answer. The raw measure beside it keeps
that factor, and therefore behaves quite differently.

### Redundancy is a count, and it does see size

Redundancy `ρ` is the raw measure that representativeness normalises, and they differ by exactly the
subcommunity's weight: `ρ̄ = w .* ρ` at the subcommunity level. Where representativeness asks what
*proportion* of the metacommunity a subcommunity stands for, redundancy asks how many times over the
things in it are represented elsewhere - a count, per individual, and so it does see size. The same
two examples make the difference concrete.

```@repl framework
raw_sub_rho(Metacommunity(twothirds), 1)[!, :diversity]
```

Every type sits in two of the three subcommunities, so lose any one and everything in it still
survives elsewhere: redundancy exactly 2. Representativeness was `2/3`, and `2 * 1/3 = 2/3`.

```@repl framework
raw_sub_rho(Metacommunity(allshared), 1)[!, :diversity]
```

**This is the case to look at.** All three subcommunities were equally and perfectly representative;
their redundancies are 6, 3 and 2 - the reciprocals of their weights. The *smallest* subcommunity is
the most redundant, because the other five sixths of the metacommunity duplicate everything in it,
while the largest holds half of everything itself and so is duplicated only twice over. Nothing about
the composition changed between the two measures; dividing by `w` is the whole of the difference.

Note: Which one you want is a real choice rather than a default. Ask "how much of the metacommunity does
this site speak for?" and you want `ρ̄`, where a small site is not penalised for being small. Ask "if I
lose this site, how much of what was in it survives?" and you want `ρ`, where being small is precisely
what makes it expendable.

### The effective number of *distinct* subcommunities is normalised in the same way

The beta measures are the reciprocals of the rho ones, and they inherit the split exactly.
[`NormalisedBeta`](@ref) `β̄` is the normalised one - the same `w` divided out - and asks how many
*distinct* subcommunities like this one the metacommunity would amount to.

```@repl framework
norm_sub_beta(Metacommunity(twothirds), 1)[!, :diversity]
norm_meta_beta(Metacommunity(twothirds), 1)[1, :diversity]
```

All three subcommunities in the matrix look much like this one, but every type is shared between two
of the three, so neither of the others is *distinct* from it - and they are worth only `1.5` distinct
subcommunities between them. Each subcommunity on its own arrives at that same number, which is why
the subcommunity-level measure is called an *estimate of* the effective number of distinct
subcommunities.

```@repl framework
norm_sub_beta(Metacommunity(allshared), 1)[!, :diversity]
norm_meta_beta(Metacommunity(allshared), 1)[1, :diversity]
```

Perfectly mixed, the answer is 1: there is effectively one **distinct** subcommunity, however
many columns the matrix has and however unequal they are. Size cannot reach this measure either - which is exactly
what lets it survive shattering, below.

### Distinctiveness is raw, and sees size the same way

[`RawBeta`](@ref) `β` is redundancy's reciprocal and the raw member of the pair, keeping the factor of
`w`. It asks what share of its own types a subcommunity holds - how much of it is its own.

```@repl framework
raw_sub_beta(Metacommunity(twothirds), 1)[!, :diversity]
raw_sub_beta(Metacommunity(allshared), 1)[!, :diversity]
```

Each type is split between two subcommunities, so each holds half of it: `0.5`, the reciprocal of the
redundancy of 2. And in the mixed case distinctiveness is `w` exactly - a subcommunity that is half
the metacommunity holds half of each of its types, and one that is a sixth holds a sixth.

Note: **Beta runs the opposite way from the other three measures**, being a reciprocal, so read it with
care: 1 is the *maximum* distinctiveness, reached only when nothing outside the subcommunity resembles
anything within it, and small values mean heavily shared. Note: That ceiling holds when types are wholly
distinct; a general similarity matrix can break it, as the section on similarity below shows.

### Gamma is per individual, so it does not track alpha

This is the distinction that makes subcommunity gamma a new measure rather than a rescaling. Here
two subcommunities draw on types that are equally abundant in the metacommunity, but the first has
twice as many of them:

```@repl framework
sizes = [1 0; 1 0; 1 0; 1 0; 0 1; 0 1] ./ 6
norm_sub_alpha(Metacommunity(sizes), 1)[!, :diversity]
sub_gamma(Metacommunity(sizes), 1)[!, :diversity]
```

The first subcommunity is twice as diverse in isolation - but each of its individuals contributes
exactly as much to the metacommunity as each of the second's, because all six types are equally
abundant overall. Rarity, not richness, is what changes the contribution:

```@repl framework
rarity = [1 0 0; 1 0 0; 0 1 2; 0 1 2] ./ 8
norm_sub_alpha(Metacommunity(rarity), 1)[!, :diversity][1:2]
sub_gamma(Metacommunity(rarity), 1)[!, :diversity][1:2]
```

Now the first two subcommunities are the same size and equally diverse in isolation, but the first
one's types are three times rarer in the metacommunity - and its contribution is three times greater.

### Hidden diversity

The case the framework was built to expose. A subcommunity holding a single very rare type is
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

The first subcommunity has `ᾱ = 1` - one type, no diversity at all - but `γ ≈ 10^9`. Adding nine more
equally rare types to the second subcommunity raises its alpha tenfold and leaves its gamma
untouched, because each individual is still just as rare; its distinctiveness falls, because a
ten-type subcommunity looks a little more like the metacommunity than a one-type one does. Alpha and
beta alone would have told you the first site was the least interesting in the study.

## Properties

**Invariance under shattering.** If you draw a boundary through a subcommunity that is internally
well mixed, you have not created a new subcommunity, and the answer should not change. The normalised
metacommunity measures - `Ā`, `R̄`, `B̄` and `G` - are invariant under such shattering:

```@repl framework
whole = [2 1; 1 3; 0 1] ./ 8
split = [2 0.5 0.5; 1 1.5 1.5; 0 0.5 0.5] ./ 8
norm_meta_beta(Metacommunity(whole), 1)[!, :diversity]
norm_meta_beta(Metacommunity(split), 1)[!, :diversity]
```

Note: The raw measures `A`, `R` and `B` are deliberately *not* invariant, and should not be: `R` is the
average redundancy of the subcommunities, and cutting one in two genuinely does create redundancy.

```@repl framework
raw_meta_rho(Metacommunity(whole), 1)[!, :diversity]
raw_meta_rho(Metacommunity(split), 1)[!, :diversity]
```

**Conditional independence.** A subcommunity's measures depend only on its own composition and on the
metacommunity as a whole - never on how the *rest* of the metacommunity happens to be divided up.
Without this, the apparent importance of one site would shift when someone redrew a boundary
elsewhere, and comparing subcommunities would be meaningless. It is what licenses the direct
comparisons above.

## Special cases

| when | what happens |
|---|---|
| `Z = I` (`UniqueTypes`) | the **naive-type** case: types wholly distinct, and the measures reduce to Hill numbers |
| no shared types between subcommunities | the **naive-community** case: `B = 1`, `B̄` is the effective number of *distinct* subcommunities, and `G = A` |
| every subcommunity has the metacommunity's composition | **well-mixed**: `R̄ = B̄ = 1`, and `Ā = G` |
| one subcommunity | every beta measure is 1, and `α = ᾱ = γ` |

The naive-type case recovers Hill numbers exactly:

```@repl framework
using Diversity.Hill
naive = [0.5, 0.3, 0.2]
meta_gamma(Metacommunity(naive), [0, 1, 2])[!, :diversity]
hillnumber(naive, [0, 1, 2])[!, :diversity]
```

### Similarity can break the bounds you expect

Read only the naive-type case and you will absorb some inequalities that do not hold in general -
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

Representativeness above 1, distinctiveness below it. Nor is this a pathological matrix - it
satisfies the triangle inequality - so no weaker assumption rescues those bounds. If your `Z` is not
the identity, check the behaviour you are relying on rather than assuming it.

```@contents
```

```@index
```
