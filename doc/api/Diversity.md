# Diversity


## Functions [Exported]

---

<a id="function__community.1" class="lexicon_definition"></a>
#### Diversity.Community [¶](#function__community.1)
### Community type, representing a single community


*source:*
[Diversity/src/Collection.jl:206](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L206)

---

<a id="function__ecosystem.1" class="lexicon_definition"></a>
#### Diversity.Ecosystem [¶](#function__ecosystem.1)
### Ecosystem type, representing an ecosystem of multiple subcommunities


*source:*
[Diversity/src/Collection.jl:198](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L198)

## Methods [Exported]

---

<a id="method__diversity.1" class="lexicon_definition"></a>
#### diversity{S<:AbstractFloat, T<:Diversity.Similarity}(kind::Tuple{Function, Set{Symbol}},  proportions::Array{S<:AbstractFloat, 2},  qs,  sim::T<:Diversity.Similarity) [¶](#method__diversity.1)
### Calculates subcommunity and supercommunity diversities

Calculates any diversity of a series of columns representing
independent subcommunity counts, for a series of orders, repesented as
a vector of qs, with similarity matrix Z, by default the (naïve)
identity matrix.

#### Arguments:
- Tuple containing:
  * `measure`: the diversity function to be used - one of Diversity.α, Diversity.ᾱ,
Diversity.ρ, Diversity.ρ̄, Diversity.γ or Diversity.γ̄
  * a Set of symbols, containing some or all of :super, :sub and :weights,
detailing what should be returned
- `proportions`:population proportions
- `qs`: single number or vector of values of parameter q
- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:
Some or all (as tuple) of:  

- vector of supercommunity diversities representing values of q  

- array of diversities, first dimension representing subcommunities, and
  last representing values of q  

- multidimensional array with dimensions matiching shape of proportions,
  with extra dimension for values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:30](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L30)

---

<a id="method__diversity.2" class="lexicon_definition"></a>
#### diversity{S<:AbstractFloat}(kind::Tuple{Function, Set{Symbol}},  proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__diversity.2)
### Calculates subcommunity and supercommunity diversities

Calculates any diversity of a series of columns representing
independent subcommunity counts, for a series of orders, repesented as
a vector of qs, with similarity matrix Z, by default the (naïve)
identity matrix.

#### Arguments:
- Tuple containing:
  * `measure`: the diversity function to be used - one of Diversity.α, Diversity.ᾱ,
Diversity.ρ, Diversity.ρ̄, Diversity.γ or Diversity.γ̄
  * a Set of symbols, containing some or all of :super, :sub and :weights,
detailing what should be returned
- `proportions`:population proportions
- `qs`: single number or vector of values of parameter q
- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:
Some or all (as tuple) of:  

- vector of supercommunity diversities representing values of q  

- array of diversities, first dimension representing subcommunities, and
  last representing values of q  

- multidimensional array with dimensions matiching shape of proportions,
  with extra dimension for values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:30](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L30)

---

<a id="method__diversity.3" class="lexicon_definition"></a>
#### diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__diversity.3)
### Calculates subcommunity and supercommunity diversities

Calculates any diversity of a series of columns representing
independent subcommunity counts, for a series of orders, repesented as
a vector of qs, with similarity matrix Z, by default the (naïve)
identity matrix.

#### Arguments:
- `measure`: the diversity function to be used - one of Dα, Dᾱ, Dρ, Dϵ
          (or Dρ̄), Dγ or Dγ̄

- `proportions`:population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

- `returnsupercommunity`: boolean describing whether to return the
                         supercommunity diversity

- `returnsubcommunity`: boolean describing whether to return the
                       subcommunity diversities 

- `returnweights`: boolean describing whether to return subcommunity weights

#### Returns:
Some or all (as tuple) of:  

- vector of supercommunity diversities representing values of q  

- array of diversities, first dimension representing subcommunities, and
  last representing values of q  

- multidimensional array with dimensions matiching shape of proportions,
  with extra dimension for values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:127](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L127)

---

<a id="method__diversity.4" class="lexicon_definition"></a>
#### diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  sim,  returnsupercommunity::Bool) [¶](#method__diversity.4)
### Calculates subcommunity and supercommunity diversities

Calculates any diversity of a series of columns representing
independent subcommunity counts, for a series of orders, repesented as
a vector of qs, with similarity matrix Z, by default the (naïve)
identity matrix.

#### Arguments:
- `measure`: the diversity function to be used - one of Dα, Dᾱ, Dρ, Dϵ
          (or Dρ̄), Dγ or Dγ̄

- `proportions`:population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

- `returnsupercommunity`: boolean describing whether to return the
                         supercommunity diversity

- `returnsubcommunity`: boolean describing whether to return the
                       subcommunity diversities 

- `returnweights`: boolean describing whether to return subcommunity weights

#### Returns:
Some or all (as tuple) of:  

- vector of supercommunity diversities representing values of q  

- array of diversities, first dimension representing subcommunities, and
  last representing values of q  

- multidimensional array with dimensions matiching shape of proportions,
  with extra dimension for values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:127](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L127)

---

<a id="method__diversity.5" class="lexicon_definition"></a>
#### diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  sim,  returnsupercommunity::Bool,  returnsubcommunity::Bool) [¶](#method__diversity.5)
### Calculates subcommunity and supercommunity diversities

Calculates any diversity of a series of columns representing
independent subcommunity counts, for a series of orders, repesented as
a vector of qs, with similarity matrix Z, by default the (naïve)
identity matrix.

#### Arguments:
- `measure`: the diversity function to be used - one of Dα, Dᾱ, Dρ, Dϵ
          (or Dρ̄), Dγ or Dγ̄

- `proportions`:population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

- `returnsupercommunity`: boolean describing whether to return the
                         supercommunity diversity

- `returnsubcommunity`: boolean describing whether to return the
                       subcommunity diversities 

- `returnweights`: boolean describing whether to return subcommunity weights

#### Returns:
Some or all (as tuple) of:  

- vector of supercommunity diversities representing values of q  

- array of diversities, first dimension representing subcommunities, and
  last representing values of q  

- multidimensional array with dimensions matiching shape of proportions,
  with extra dimension for values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:127](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L127)

---

<a id="method__diversity.6" class="lexicon_definition"></a>
#### diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  sim,  returnsupercommunity::Bool,  returnsubcommunity::Bool,  returnweights::Bool) [¶](#method__diversity.6)
### Calculates subcommunity and supercommunity diversities

Calculates any diversity of a series of columns representing
independent subcommunity counts, for a series of orders, repesented as
a vector of qs, with similarity matrix Z, by default the (naïve)
identity matrix.

#### Arguments:
- `measure`: the diversity function to be used - one of Dα, Dᾱ, Dρ, Dϵ
          (or Dρ̄), Dγ or Dγ̄

- `proportions`:population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

- `returnsupercommunity`: boolean describing whether to return the
                         supercommunity diversity

- `returnsubcommunity`: boolean describing whether to return the
                       subcommunity diversities 

- `returnweights`: boolean describing whether to return subcommunity weights

#### Returns:
Some or all (as tuple) of:  

- vector of supercommunity diversities representing values of q  

- array of diversities, first dimension representing subcommunities, and
  last representing values of q  

- multidimensional array with dimensions matiching shape of proportions,
  with extra dimension for values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:127](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L127)

---

<a id="method__qdz.1" class="lexicon_definition"></a>
#### qDZ{S<:AbstractFloat, T<:Diversity.Similarity}(proportions::Array{S<:AbstractFloat, 1},  qs,  sim::T<:Diversity.Similarity) [¶](#method__qdz.1)
### Calculates Leinster-Cobbold / similarity-sensitive diversity

Calculates Leinster-Cobbold general diversity of >= 1 order(s) *qs* of
a population with given relative *proportions*, and similarity matrix
*Z*.

#### Arguments:
- `proportions`: relative proportions of different individuals /
types in a population or series of populations
- `qs`: single number or vector of orders of diversity measurement
- `Z`: similarity matrix

#### Returns:
- Diversity of order qs (single number or vector of diversities)


*source:*
[Diversity/src/EffectiveNumbers.jl:113](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/EffectiveNumbers.jl#L113)

---

<a id="method__qdz.2" class="lexicon_definition"></a>
#### qDZ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 1},  qs) [¶](#method__qdz.2)
### Calculates Leinster-Cobbold / similarity-sensitive diversity

Calculates Leinster-Cobbold general diversity of >= 1 order(s) *qs* of
a population with given relative *proportions*, and similarity matrix
*Z*.

#### Arguments:
- `proportions`: relative proportions of different individuals /
types in a population or series of populations
- `qs`: single number or vector of orders of diversity measurement
- `Z`: similarity matrix

#### Returns:
- Diversity of order qs (single number or vector of diversities)


*source:*
[Diversity/src/EffectiveNumbers.jl:113](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/EffectiveNumbers.jl#L113)

---

<a id="method__qd.1" class="lexicon_definition"></a>
#### qD{S<:Real}(proportions::Array{S<:Real, 1},  qs) [¶](#method__qd.1)
### Calculates Hill / naive-similarity diversity

Calculates Hill number or naive diversity of order(s) *qs* of a
population with given relative proportions.

#### Arguments:
- `proportions`: relative proportions of different individuals /
types in population or series of populations
- `qs`: single number or vector of orders of diversity measurement

#### Returns:
- Diversity of order qs (single number or vector of diversities)

*source:*
[Diversity/src/EffectiveNumbers.jl:85](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/EffectiveNumbers.jl#L85)

---

<a id="method__subcommunityalphabar.1" class="lexicon_definition"></a>
#### subcommunityalphabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__subcommunityalphabar.1)
### Normalised similarity-sensitive subcommunity alpha diversity)

Calculates (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:222](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L222)

---

<a id="method__subcommunityalphabar.2" class="lexicon_definition"></a>
#### subcommunityalphabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__subcommunityalphabar.2)
### Normalised similarity-sensitive subcommunity alpha diversity)

Calculates (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:222](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L222)

---

<a id="method__subcommunityalpha.1" class="lexicon_definition"></a>
#### subcommunityalpha{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__subcommunityalpha.1)
### Raw similarity-sensitive subcommunity alpha diversity / naive-community diversity

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:198](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L198)

---

<a id="method__subcommunityalpha.2" class="lexicon_definition"></a>
#### subcommunityalpha{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__subcommunityalpha.2)
### Raw similarity-sensitive subcommunity alpha diversity / naive-community diversity

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:198](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L198)

---

<a id="method__subcommunitybetabar.1" class="lexicon_definition"></a>
#### subcommunitybetabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__subcommunitybetabar.1)
### Normalised similarity-sensitive subcommunity beta diversity

Calculates normalised beta diversities or the effective number of
distinct subcommunities perceived by a series of subcommunities
represented by columns of independent subcommunity counts, represented
as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:371](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L371)

---

<a id="method__subcommunitybetabar.2" class="lexicon_definition"></a>
#### subcommunitybetabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__subcommunitybetabar.2)
### Normalised similarity-sensitive subcommunity beta diversity

Calculates normalised beta diversities or the effective number of
distinct subcommunities perceived by a series of subcommunities
represented by columns of independent subcommunity counts, represented
as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:371](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L371)

---

<a id="method__subcommunitybeta.1" class="lexicon_definition"></a>
#### subcommunitybeta{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__subcommunitybeta.1)
### Raw similarity-sensitive subcommunity beta diversity / distinctiveness / concentration

Calculates the raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:316](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L316)

---

<a id="method__subcommunitybeta.2" class="lexicon_definition"></a>
#### subcommunitybeta{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__subcommunitybeta.2)
### Raw similarity-sensitive subcommunity beta diversity / distinctiveness / concentration

Calculates the raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:316](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L316)

---

<a id="method__subcommunitygammabar.1" class="lexicon_definition"></a>
#### subcommunitygammabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__subcommunitygammabar.1)
### Normalised similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:517](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L517)

---

<a id="method__subcommunitygammabar.2" class="lexicon_definition"></a>
#### subcommunitygammabar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__subcommunitygammabar.2)
### Normalised similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:517](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L517)

---

<a id="method__subcommunitygamma.1" class="lexicon_definition"></a>
#### subcommunitygamma{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__subcommunitygamma.1)
### Raw similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:494](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L494)

---

<a id="method__subcommunitygamma.2" class="lexicon_definition"></a>
#### subcommunitygamma{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__subcommunitygamma.2)
### Raw similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:494](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L494)

---

<a id="method__subcommunityrhobar.1" class="lexicon_definition"></a>
#### subcommunityrhobar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__subcommunityrhobar.1)
### Normalised similarity-sensitive subcommunity representativeness

Calculates the representativeness of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs. Representativeness
reflects what proportion of the supercommunity each subcommunity is
representative of on average, so if each subcommunity contains 1/xth
of the species, then the average representativeness of the
subcommunities is 1/x.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of representativenesses, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:345](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L345)

---

<a id="method__subcommunityrhobar.2" class="lexicon_definition"></a>
#### subcommunityrhobar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__subcommunityrhobar.2)
### Normalised similarity-sensitive subcommunity representativeness

Calculates the representativeness of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs. Representativeness
reflects what proportion of the supercommunity each subcommunity is
representative of on average, so if each subcommunity contains 1/xth
of the species, then the average representativeness of the
subcommunities is 1/x.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of representativenesses, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:345](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L345)

---

<a id="method__subcommunityrho.1" class="lexicon_definition"></a>
#### subcommunityrho{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__subcommunityrho.1)
### Raw similarity-sensitive subcommunity redundancy

Calculates redundancy of a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of redundancies, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:290](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L290)

---

<a id="method__subcommunityrho.2" class="lexicon_definition"></a>
#### subcommunityrho{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__subcommunityrho.2)
### Raw similarity-sensitive subcommunity redundancy

Calculates redundancy of a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of redundancies, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:290](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L290)

---

<a id="method__supercommunityabar.1" class="lexicon_definition"></a>
#### supercommunityAbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__supercommunityabar.1)
### Normalised similarity-sensitive supercommunity alpha diversity

Calculates average (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:267](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L267)

---

<a id="method__supercommunityabar.2" class="lexicon_definition"></a>
#### supercommunityAbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__supercommunityabar.2)
### Normalised similarity-sensitive supercommunity alpha diversity

Calculates average (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:267](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L267)

---

<a id="method__supercommunitya.1" class="lexicon_definition"></a>
#### supercommunityA{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__supercommunitya.1)
### Raw similarity-sensitive supercommunity alpha diversity / naive-community diversity

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector
of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:245](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L245)

---

<a id="method__supercommunitya.2" class="lexicon_definition"></a>
#### supercommunityA{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__supercommunitya.2)
### Raw similarity-sensitive supercommunity alpha diversity / naive-community diversity

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector
of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:245](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L245)

---

<a id="method__supercommunitybbar.1" class="lexicon_definition"></a>
#### supercommunityBbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__supercommunitybbar.1)
### Normalised similarity-sensitive supercommunity beta diversity / effective number of communities

Calculates average normalised beta diversity or the effective number
of distinct subcommunities present in a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:471](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L471)

---

<a id="method__supercommunitybbar.2" class="lexicon_definition"></a>
#### supercommunityBbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__supercommunitybbar.2)
### Normalised similarity-sensitive supercommunity beta diversity / effective number of communities

Calculates average normalised beta diversity or the effective number
of distinct subcommunities present in a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:471](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L471)

---

<a id="method__supercommunityb.1" class="lexicon_definition"></a>
#### supercommunityB{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__supercommunityb.1)
### Raw similarity-sensitive supercommunity beta diversity / distinctiveness / concentration

Calculates average raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:418](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L418)

---

<a id="method__supercommunityb.2" class="lexicon_definition"></a>
#### supercommunityB{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__supercommunityb.2)
### Raw similarity-sensitive supercommunity beta diversity / distinctiveness / concentration

Calculates average raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:418](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L418)

---

<a id="method__supercommunitygbar.1" class="lexicon_definition"></a>
#### supercommunityGbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__supercommunitygbar.1)
### Normalised similarity-sensitive supercommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:561](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L561)

---

<a id="method__supercommunitygbar.2" class="lexicon_definition"></a>
#### supercommunityGbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__supercommunitygbar.2)
### Normalised similarity-sensitive supercommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:561](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L561)

---

<a id="method__supercommunityg.1" class="lexicon_definition"></a>
#### supercommunityG{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__supercommunityg.1)
### Raw similarity-sensitive supercommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:539](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L539)

---

<a id="method__supercommunityg.2" class="lexicon_definition"></a>
#### supercommunityG{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__supercommunityg.2)
### Raw similarity-sensitive supercommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:539](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L539)

---

<a id="method__supercommunityrbar.1" class="lexicon_definition"></a>
#### supercommunityRbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__supercommunityrbar.1)
### Normalised similarity-sensitive supercommunity representativeness

Calculates average representativeness of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs. Representativeness
reflects what proportion of the supercommunity each subcommunity is
representative of on average, so if each subcommunity contains 1/xth
of the species, then the average representativeness of the
subcommunities is 1/x.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of representativenesses representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:446](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L446)

---

<a id="method__supercommunityrbar.2" class="lexicon_definition"></a>
#### supercommunityRbar{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__supercommunityrbar.2)
### Normalised similarity-sensitive supercommunity representativeness

Calculates average representativeness of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs. Representativeness
reflects what proportion of the supercommunity each subcommunity is
representative of on average, so if each subcommunity contains 1/xth
of the species, then the average representativeness of the
subcommunities is 1/x.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of representativenesses representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:446](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L446)

---

<a id="method__supercommunityr.1" class="lexicon_definition"></a>
#### supercommunityR{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__supercommunityr.1)
### Raw similarity-sensitive supercommunity redundancy

Calculates average redundancy of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of redundancies representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:393](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L393)

---

<a id="method__supercommunityr.2" class="lexicon_definition"></a>
#### supercommunityR{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  sim) [¶](#method__supercommunityr.2)
### Raw similarity-sensitive supercommunity redundancy

Calculates average redundancy of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of redundancies representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:393](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/GeneralisedDiversities.jl#L393)

## Types [Exported]

---

<a id="type__collection.1" class="lexicon_definition"></a>
#### Diversity.Collection{S<:Diversity.Similarity, P<:Diversity.Partition, FP<:AbstractFloat} [¶](#type__collection.1)
### Collection type, representing a collection of one or more subcommunities

Type representing a single community or collection of communities. It
contains a collection of individuals which *may* be further
partitioned into smaller groups. For instance this may be an
ecosystem, which consists of a series of subcommunities.

The type stores relative abundances of different types, e.g. species,
and also allows for similarity between individuals.

#### Parameterisation:

**Collection{S, P, FP}**

- `S` is the similarity type, e.g. Species, a subtype of Similarity.

- `P` is the partition type, e.g. Subcommunity, a subtype of Partition.

- `FP` is the kind of number storage, a subtype of AbstractFloat.

#### Members:

- `abundances` An array of relative abundances. The first dimension
               represents the species, and further dimensions
               represent the structure of collection.

- `z` A two-dimensional matrix representing similarity between
      individuals of the base type, S. By default this will be the
      identity matrix.


*source:*
[Diversity/src/Collection.jl:147](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L147)

---

<a id="type__generalsimilarity.1" class="lexicon_definition"></a>
#### Diversity.GeneralSimilarity{S<:AbstractFloat} [¶](#type__generalsimilarity.1)
### A general matrix-based Similarity subtype

This subtype of Similarity simply holds a matrix with similarities
between individuals.

#### Members:

- `matrix` A two-dimensional matrix representing similarity between
    individuals. By default this will be the identity matrix,
    but will require the number of species to be instantiated.


*source:*
[Diversity/src/Collection.jl:48](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L48)

---

<a id="type__onecommunity.1" class="lexicon_definition"></a>
#### Diversity.Onecommunity [¶](#type__onecommunity.1)
### Partition type allowing only one subcommunity


*source:*
[Diversity/src/Collection.jl:108](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L108)

---

<a id="type__subcommunity.1" class="lexicon_definition"></a>
#### Diversity.Subcommunity [¶](#type__subcommunity.1)
### Partition type with multiple subccomunities


*source:*
[Diversity/src/Collection.jl:103](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L103)

---

<a id="type__taxonomy.1" class="lexicon_definition"></a>
#### Diversity.Taxonomy [¶](#type__taxonomy.1)
### A subtype of Similarity with similarity between related taxa

This subtype of Similarity allows taxonomic similarity matrices


*source:*
[Diversity/src/Collection.jl:32](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L32)

---

<a id="type__unique.1" class="lexicon_definition"></a>
#### Diversity.Unique [¶](#type__unique.1)
### A subtype of Similarity where all individuals are completely distinct

This type is the simplest Similarity subtype, which identifies all
individuals as unique and completely distinct from each other.


*source:*
[Diversity/src/Collection.jl:16](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L16)

## Typealiass [Exported]

---

<a id="typealias__species.1" class="lexicon_definition"></a>
#### Species [¶](#typealias__species.1)
### A subtype of Similarity where all species are completely distinct

This type is the simplest Similarity subtype, which identifies all
species as unique and completely distinct from each other.


*source:*
[Diversity/src/Collection.jl:25](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L25)


## Methods [Internal]

---

<a id="method__call.1" class="lexicon_definition"></a>
#### call{S<:AbstractFloat}(::Type{Diversity.GeneralSimilarity{S<:AbstractFloat}},  z::Array{S<:AbstractFloat, 2}) [¶](#method__call.1)
### Constructor for GeneralSimilarity

Creates an instance of the GeneralSimilarity class, with an arbitrary similarity matrix.

#### Arguments:
- `z`: similarity matrix


*source:*
[Diversity/src/Collection.jl:65](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L65)

---

<a id="method__contributions.1" class="lexicon_definition"></a>
#### contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__contributions.1)
### Calculate diversity contributions from subcommunities

Calculates proportions that subcommunities each contribute to
supercommunity diversity per subcommunity (perindividual = false), or
per individual (perindividual = true) - in the latter case scaled
so that the total # of individuals is 1, since we only have
relative abundances.

#### Arguments:
- `measure`: diversity measure to use
- `proportions`: population proportions
- `qs`: single number or vector of values of parameter q
- `perindividual`: do we measure per individual in population (true)
                   or per subcommunity (false)
- `Z`: similarity matrix
- `returnsupercommunity`: boolean describing whether to return the
                    supercommunity diversity
- `returnsubcommunity`: boolean describing whether to return the
                    subcommunity diversities
- `returnweights`: boolean describing whether to return subcommunity weights

#### Returns:
- contributions of subcommunities to supercommunity diversity (of type measure)
- and none, some or all (in a tuple) of:
  - vector of supercommunity diversities representing values of q
  - array of diversities, first dimension representing subcommunities, and
    last representing values of q
  - vector of subcommunity weights


*source:*
[Diversity/src/CommunityContributions.jl:31](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/CommunityContributions.jl#L31)

---

<a id="method__contributions.2" class="lexicon_definition"></a>
#### contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  perindividual::Bool) [¶](#method__contributions.2)
### Calculate diversity contributions from subcommunities

Calculates proportions that subcommunities each contribute to
supercommunity diversity per subcommunity (perindividual = false), or
per individual (perindividual = true) - in the latter case scaled
so that the total # of individuals is 1, since we only have
relative abundances.

#### Arguments:
- `measure`: diversity measure to use
- `proportions`: population proportions
- `qs`: single number or vector of values of parameter q
- `perindividual`: do we measure per individual in population (true)
                   or per subcommunity (false)
- `Z`: similarity matrix
- `returnsupercommunity`: boolean describing whether to return the
                    supercommunity diversity
- `returnsubcommunity`: boolean describing whether to return the
                    subcommunity diversities
- `returnweights`: boolean describing whether to return subcommunity weights

#### Returns:
- contributions of subcommunities to supercommunity diversity (of type measure)
- and none, some or all (in a tuple) of:
  - vector of supercommunity diversities representing values of q
  - array of diversities, first dimension representing subcommunities, and
    last representing values of q
  - vector of subcommunity weights


*source:*
[Diversity/src/CommunityContributions.jl:31](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/CommunityContributions.jl#L31)

---

<a id="method__contributions.3" class="lexicon_definition"></a>
#### contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  perindividual::Bool,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__contributions.3)
### Calculate diversity contributions from subcommunities

Calculates proportions that subcommunities each contribute to
supercommunity diversity per subcommunity (perindividual = false), or
per individual (perindividual = true) - in the latter case scaled
so that the total # of individuals is 1, since we only have
relative abundances.

#### Arguments:
- `measure`: diversity measure to use
- `proportions`: population proportions
- `qs`: single number or vector of values of parameter q
- `perindividual`: do we measure per individual in population (true)
                   or per subcommunity (false)
- `Z`: similarity matrix
- `returnsupercommunity`: boolean describing whether to return the
                    supercommunity diversity
- `returnsubcommunity`: boolean describing whether to return the
                    subcommunity diversities
- `returnweights`: boolean describing whether to return subcommunity weights

#### Returns:
- contributions of subcommunities to supercommunity diversity (of type measure)
- and none, some or all (in a tuple) of:
  - vector of supercommunity diversities representing values of q
  - array of diversities, first dimension representing subcommunities, and
    last representing values of q
  - vector of subcommunity weights


*source:*
[Diversity/src/CommunityContributions.jl:31](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/CommunityContributions.jl#L31)

---

<a id="method__contributions.4" class="lexicon_definition"></a>
#### contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  perindividual::Bool,  Z::Array{S<:AbstractFloat, 2},  returnsupercommunity::Bool) [¶](#method__contributions.4)
### Calculate diversity contributions from subcommunities

Calculates proportions that subcommunities each contribute to
supercommunity diversity per subcommunity (perindividual = false), or
per individual (perindividual = true) - in the latter case scaled
so that the total # of individuals is 1, since we only have
relative abundances.

#### Arguments:
- `measure`: diversity measure to use
- `proportions`: population proportions
- `qs`: single number or vector of values of parameter q
- `perindividual`: do we measure per individual in population (true)
                   or per subcommunity (false)
- `Z`: similarity matrix
- `returnsupercommunity`: boolean describing whether to return the
                    supercommunity diversity
- `returnsubcommunity`: boolean describing whether to return the
                    subcommunity diversities
- `returnweights`: boolean describing whether to return subcommunity weights

#### Returns:
- contributions of subcommunities to supercommunity diversity (of type measure)
- and none, some or all (in a tuple) of:
  - vector of supercommunity diversities representing values of q
  - array of diversities, first dimension representing subcommunities, and
    last representing values of q
  - vector of subcommunity weights


*source:*
[Diversity/src/CommunityContributions.jl:31](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/CommunityContributions.jl#L31)

---

<a id="method__contributions.5" class="lexicon_definition"></a>
#### contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  perindividual::Bool,  Z::Array{S<:AbstractFloat, 2},  returnsupercommunity::Bool,  returnsubcommunity::Bool) [¶](#method__contributions.5)
### Calculate diversity contributions from subcommunities

Calculates proportions that subcommunities each contribute to
supercommunity diversity per subcommunity (perindividual = false), or
per individual (perindividual = true) - in the latter case scaled
so that the total # of individuals is 1, since we only have
relative abundances.

#### Arguments:
- `measure`: diversity measure to use
- `proportions`: population proportions
- `qs`: single number or vector of values of parameter q
- `perindividual`: do we measure per individual in population (true)
                   or per subcommunity (false)
- `Z`: similarity matrix
- `returnsupercommunity`: boolean describing whether to return the
                    supercommunity diversity
- `returnsubcommunity`: boolean describing whether to return the
                    subcommunity diversities
- `returnweights`: boolean describing whether to return subcommunity weights

#### Returns:
- contributions of subcommunities to supercommunity diversity (of type measure)
- and none, some or all (in a tuple) of:
  - vector of supercommunity diversities representing values of q
  - array of diversities, first dimension representing subcommunities, and
    last representing values of q
  - vector of subcommunity weights


*source:*
[Diversity/src/CommunityContributions.jl:31](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/CommunityContributions.jl#L31)

---

<a id="method__contributions.6" class="lexicon_definition"></a>
#### contributions{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  perindividual::Bool,  Z::Array{S<:AbstractFloat, 2},  returnsupercommunity::Bool,  returnsubcommunity::Bool,  returnweights::Bool) [¶](#method__contributions.6)
### Calculate diversity contributions from subcommunities

Calculates proportions that subcommunities each contribute to
supercommunity diversity per subcommunity (perindividual = false), or
per individual (perindividual = true) - in the latter case scaled
so that the total # of individuals is 1, since we only have
relative abundances.

#### Arguments:
- `measure`: diversity measure to use
- `proportions`: population proportions
- `qs`: single number or vector of values of parameter q
- `perindividual`: do we measure per individual in population (true)
                   or per subcommunity (false)
- `Z`: similarity matrix
- `returnsupercommunity`: boolean describing whether to return the
                    supercommunity diversity
- `returnsubcommunity`: boolean describing whether to return the
                    subcommunity diversities
- `returnweights`: boolean describing whether to return subcommunity weights

#### Returns:
- contributions of subcommunities to supercommunity diversity (of type measure)
- and none, some or all (in a tuple) of:
  - vector of supercommunity diversities representing values of q
  - array of diversities, first dimension representing subcommunities, and
    last representing values of q
  - vector of subcommunity weights


*source:*
[Diversity/src/CommunityContributions.jl:31](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/CommunityContributions.jl#L31)

---

<a id="method__powermean.1" class="lexicon_definition"></a>
#### powermean{S<:AbstractFloat}(values::Array{S<:AbstractFloat, 1}) [¶](#method__powermean.1)
### Calculates the weighted powermean of a series of numbers

Calculates *order*th power mean of *values*, weighted by
*weights*. By default, *weights* are equal and *order*
is 1, so this is just the arithmetic mean.

#### Arguments:
- `values`: values for which to calculate mean
- `order`: order of power mean
- `weights`: weights of elements, normalised to 1 inside function

#### Returns:
- weighted power mean(s)


*source:*
[Diversity/src/EffectiveNumbers.jl:16](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/EffectiveNumbers.jl#L16)

---

<a id="method__powermean.2" class="lexicon_definition"></a>
#### powermean{S<:AbstractFloat}(values::Array{S<:AbstractFloat, 1},  order::S<:AbstractFloat) [¶](#method__powermean.2)
### Calculates the weighted powermean of a series of numbers

Calculates *order*th power mean of *values*, weighted by
*weights*. By default, *weights* are equal and *order*
is 1, so this is just the arithmetic mean.

#### Arguments:
- `values`: values for which to calculate mean
- `order`: order of power mean
- `weights`: weights of elements, normalised to 1 inside function

#### Returns:
- weighted power mean(s)


*source:*
[Diversity/src/EffectiveNumbers.jl:16](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/EffectiveNumbers.jl#L16)

---

<a id="method__powermean.3" class="lexicon_definition"></a>
#### powermean{S<:AbstractFloat}(values::Array{S<:AbstractFloat, 1},  order::S<:AbstractFloat,  weights::Array{S<:AbstractFloat, 1}) [¶](#method__powermean.3)
### Calculates the weighted powermean of a series of numbers

Calculates *order*th power mean of *values*, weighted by
*weights*. By default, *weights* are equal and *order*
is 1, so this is just the arithmetic mean.

#### Arguments:
- `values`: values for which to calculate mean
- `order`: order of power mean
- `weights`: weights of elements, normalised to 1 inside function

#### Returns:
- weighted power mean(s)


*source:*
[Diversity/src/EffectiveNumbers.jl:16](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/EffectiveNumbers.jl#L16)

## Types [Internal]

---

<a id="type__partition.1" class="lexicon_definition"></a>
#### Diversity.Partition [¶](#type__partition.1)
### Abstract Partition supertype for all partitioning types

This type is the abstract superclass of all partitioning types.
Partition subtypes allow you to define how to partition your total
collection (e.g. an ecosystem) into smaller components (e.g.
subcommunities).


*source:*
[Diversity/src/Collection.jl:98](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L98)

---

<a id="type__similarity.1" class="lexicon_definition"></a>
#### Diversity.Similarity [¶](#type__similarity.1)
### Abstract Similarity supertype for all similarity measures

This type is the abstract superclass of all similarity types. Its
subtypes allow you to define how similarity is measured between
individuals.


*source:*
[Diversity/src/Collection.jl:8](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Collection.jl#L8)

