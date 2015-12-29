# Diversity.Jost


## Methods [Exported]

---

<a id="method__jostalpha.1" class="lexicon_definition"></a>
#### jostalpha{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__jostalpha.1)
### Calculates Jost's alpha diversity

Calculates Jost's alpha diversity of a series of columns representing
independent community counts, for a series of orders, repesented as a
vector of qs. This is just the naive-community ecosystem diversity
divided by the naive-community beta diversity.

### Arguments:
- `proportions` relative proportions of different individuals / species
                in population (vector, or matrix where columns are
                for individual sub-communities)

- `qs` single number or vector of orders of diversity measurement

### Returns:
- array of diversities, first dimension representing sub-communities, and
  last representing values of q


*source:*
[Diversity/src/Jost.jl:22](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Jost.jl#L22)

---

<a id="method__jostbeta.1" class="lexicon_definition"></a>
#### jostbeta{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__jostbeta.1)
### Calculates Jost's beta diversity

Calculates Jost's beta diversity of a series of columns representing
independent community counts, for a series of orders, repesented as a
vector of qs. This is just the naive gamma diversity divided by
Jost's alpha diversity

### Arguments:
- `proportions` relative proportions of different individuals / species
                in population (vector, or matrix where columns are
                for individual sub-communities)

- `qs` single number or vector of orders of diversity measurement

### Returns:
- array of diversities, first dimension representing sub-communities, and
  last representing values of q


*source:*
[Diversity/src/Jost.jl:48](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Jost.jl#L48)

