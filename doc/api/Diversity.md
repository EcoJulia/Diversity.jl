# Diversity


## Functions [Exported]

---

<a id="function__community.1" class="lexicon_definition"></a>
#### Diversity.Community [¶](#function__community.1)
### Community type, representing a single community


*source:*
[Diversity/src/Collection.jl:213](file:///Users/richardr/.julia/v0.4/Diversity/src/Collection.jl)

---

<a id="function__ecosystem.1" class="lexicon_definition"></a>
#### Diversity.Ecosystem [¶](#function__ecosystem.1)
### Ecosystem type, representing an ecosystem of multiple subcommunities


*source:*
[Diversity/src/Collection.jl:205](file:///Users/richardr/.julia/v0.4/Diversity/src/Collection.jl)

## Methods [Exported]

---

<a id="method__da.1" class="lexicon_definition"></a>
#### DA{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__da.1)
### Raw similarity-sensitive supercommunity alpha diversity / naive-community diversity

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector
of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:162](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__da.2" class="lexicon_definition"></a>
#### DA{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__da.2)
### Raw similarity-sensitive supercommunity alpha diversity / naive-community diversity

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector
of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:162](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__db.1" class="lexicon_definition"></a>
#### DB{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__db.1)
### Raw similarity-sensitive supercommunity beta diversity / distinctiveness / concentration

Calculates average raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:378](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__db.2" class="lexicon_definition"></a>
#### DB{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__db.2)
### Raw similarity-sensitive supercommunity beta diversity / distinctiveness / concentration

Calculates average raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:378](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__db772.1" class="lexicon_definition"></a>
#### DB̄{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__db772.1)
### Normalised similarity-sensitive supercommunity beta diversity / effective number of communities

Calculates average normalised beta diversity or the effective number
of distinct subcommunities present in a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:433](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__db772.2" class="lexicon_definition"></a>
#### DB̄{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__db772.2)
### Normalised similarity-sensitive supercommunity beta diversity / effective number of communities

Calculates average normalised beta diversity or the effective number
of distinct subcommunities present in a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:433](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__dg.1" class="lexicon_definition"></a>
#### DG{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__dg.1)
### Raw similarity-sensitive supercommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:524](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__dg.2" class="lexicon_definition"></a>
#### DG{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__dg.2)
### Raw similarity-sensitive supercommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:524](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__dr.1" class="lexicon_definition"></a>
#### DR{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__dr.1)
### Raw similarity-sensitive supercommunity redundancy

Calculates average redundancy of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of redundancies representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:353](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__dr.2" class="lexicon_definition"></a>
#### DR{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__dr.2)
### Raw similarity-sensitive supercommunity redundancy

Calculates average redundancy of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of redundancies representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:353](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__dr772.1" class="lexicon_definition"></a>
#### DR̄{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__dr772.1)
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

- `Z`: similarity matrix

#### Returns:

- vector of representativenesses representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:407](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__dr772.2" class="lexicon_definition"></a>
#### DR̄{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__dr772.2)
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

- `Z`: similarity matrix

#### Returns:

- vector of representativenesses representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:407](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d256.1" class="lexicon_definition"></a>
#### DĀ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__d256.1)
### Normalised similarity-sensitive supercommunity alpha diversity

Calculates average (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:187](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d256.2" class="lexicon_definition"></a>
#### DĀ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__d256.2)
### Normalised similarity-sensitive supercommunity alpha diversity

Calculates average (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:187](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d945.1" class="lexicon_definition"></a>
#### Dα{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__d945.1)
### Raw similarity-sensitive subcommunity alpha diversity / naive-community diversity

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:100](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d945.2" class="lexicon_definition"></a>
#### Dα{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__d945.2)
### Raw similarity-sensitive subcommunity alpha diversity / naive-community diversity

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:100](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d946.1" class="lexicon_definition"></a>
#### Dβ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__d946.1)
### Raw similarity-sensitive subcommunity beta diversity / distinctiveness / concentration

Calculates the raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:246](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d946.2" class="lexicon_definition"></a>
#### Dβ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__d946.2)
### Raw similarity-sensitive subcommunity beta diversity / distinctiveness / concentration

Calculates the raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:246](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d946772.1" class="lexicon_definition"></a>
#### Dβ̄{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__d946772.1)
### Normalised similarity-sensitive subcommunity beta diversity

Calculates normalised beta diversities or the effective number of
distinct subcommunities perceived by a series of subcommunities
represented by columns of independent subcommunity counts, represented
as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:319](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d946772.2" class="lexicon_definition"></a>
#### Dβ̄{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__d946772.2)
### Normalised similarity-sensitive subcommunity beta diversity

Calculates normalised beta diversities or the effective number of
distinct subcommunities perceived by a series of subcommunities
represented by columns of independent subcommunity counts, represented
as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:319](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d947.1" class="lexicon_definition"></a>
#### Dγ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__d947.1)
### Raw similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:459](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d947.2" class="lexicon_definition"></a>
#### Dγ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__d947.2)
### Raw similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:459](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d947772.1" class="lexicon_definition"></a>
#### Dγ̄{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__d947772.1)
### Normalised similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:492](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d947772.2" class="lexicon_definition"></a>
#### Dγ̄{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__d947772.2)
### Normalised similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:492](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d961.1" class="lexicon_definition"></a>
#### Dρ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__d961.1)
### Raw similarity-sensitive subcommunity redundancy

Calculates redundancy of a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of redundancies, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:213](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d961.2" class="lexicon_definition"></a>
#### Dρ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__d961.2)
### Raw similarity-sensitive subcommunity redundancy

Calculates redundancy of a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of redundancies, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:213](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d961772.1" class="lexicon_definition"></a>
#### Dρ̄{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__d961772.1)
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

- `Z`: similarity matrix

#### Returns:

- array of representativenesses, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:283](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d961772.2" class="lexicon_definition"></a>
#### Dρ̄{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__d961772.2)
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

- `Z`: similarity matrix

#### Returns:

- array of representativenesses, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:283](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d7712.1" class="lexicon_definition"></a>
#### DḠ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__d7712.1)
### Normalised similarity-sensitive supercommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:549](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d7712.2" class="lexicon_definition"></a>
#### DḠ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__d7712.2)
### Normalised similarity-sensitive supercommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- vector of diversities representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:549](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d8113.1" class="lexicon_definition"></a>
#### Dᾱ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__d8113.1)
### Normalised similarity-sensitive subcommunity alpha diversity)

Calculates (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:131](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__d8113.2" class="lexicon_definition"></a>
#### Dᾱ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__d8113.2)
### Normalised similarity-sensitive subcommunity alpha diversity)

Calculates (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q


*source:*
[Diversity/src/GeneralisedDiversities.jl:131](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__diversity.1" class="lexicon_definition"></a>
#### diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs) [¶](#method__diversity.1)
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

- `Z`: similarity matrix

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
[Diversity/src/GeneralisedDiversities.jl:38](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__diversity.2" class="lexicon_definition"></a>
#### diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__diversity.2)
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

- `Z`: similarity matrix

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
[Diversity/src/GeneralisedDiversities.jl:38](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__diversity.3" class="lexicon_definition"></a>
#### diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2},  returnsupercommunity::Bool) [¶](#method__diversity.3)
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

- `Z`: similarity matrix

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
[Diversity/src/GeneralisedDiversities.jl:38](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__diversity.4" class="lexicon_definition"></a>
#### diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2},  returnsupercommunity::Bool,  returnsubcommunity::Bool) [¶](#method__diversity.4)
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

- `Z`: similarity matrix

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
[Diversity/src/GeneralisedDiversities.jl:38](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__diversity.5" class="lexicon_definition"></a>
#### diversity{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  qs,  Z::Array{S<:AbstractFloat, 2},  returnsupercommunity::Bool,  returnsubcommunity::Bool,  returnweights::Bool) [¶](#method__diversity.5)
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

- `Z`: similarity matrix

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
[Diversity/src/GeneralisedDiversities.jl:38](file:///Users/richardr/.julia/v0.4/Diversity/src/GeneralisedDiversities.jl)

---

<a id="method__qdz.1" class="lexicon_definition"></a>
#### qDZ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 1},  qs) [¶](#method__qdz.1)
### Calculates Leinster-Cobbold / similarity-sensitive diversity

Calculates Leinster-Cobbold general diversity of >= 1 order(s) *qs* of
a population with given relative *proportions*, and similarity matrix
*Z*.

#### Arguments:
- `proportions`: relative proportions of different individuals /
               species in a population or series of populations
- `qs`: single number or vector of orders of diversity measurement
- `Z`: similarity matrix

#### Returns:
- Diversity of order qs (single number or vector of diversities)

*source:*
[Diversity/src/EffectiveNumbers.jl:112](file:///Users/richardr/.julia/v0.4/Diversity/src/EffectiveNumbers.jl)

---

<a id="method__qdz.2" class="lexicon_definition"></a>
#### qDZ{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 1},  qs,  Z::Array{S<:AbstractFloat, 2}) [¶](#method__qdz.2)
### Calculates Leinster-Cobbold / similarity-sensitive diversity

Calculates Leinster-Cobbold general diversity of >= 1 order(s) *qs* of
a population with given relative *proportions*, and similarity matrix
*Z*.

#### Arguments:
- `proportions`: relative proportions of different individuals /
               species in a population or series of populations
- `qs`: single number or vector of orders of diversity measurement
- `Z`: similarity matrix

#### Returns:
- Diversity of order qs (single number or vector of diversities)

*source:*
[Diversity/src/EffectiveNumbers.jl:112](file:///Users/richardr/.julia/v0.4/Diversity/src/EffectiveNumbers.jl)

---

<a id="method__qd.1" class="lexicon_definition"></a>
#### qD{S<:Number}(proportions::Array{S<:Number, 1},  qs) [¶](#method__qd.1)
### Calculates Hill / naive-similarity diversity

Calculates Hill number or naive diversity of order(s) *qs* of a
population with given relative proportions.

#### Arguments:
- `proportions`: relative proportions of different individuals /
               species in population or series of populations
- `qs`: single number or vector of orders of diversity measurement

#### Returns:
- Diversity of order qs (single number or vector of diversities)

*source:*
[Diversity/src/EffectiveNumbers.jl:85](file:///Users/richardr/.julia/v0.4/Diversity/src/EffectiveNumbers.jl)

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

- `Z` A two-dimensional matrix representing similarity between
      individuals of the base type, S. By default this will be the
      identity matrix.


*source:*
[Diversity/src/Collection.jl:154](file:///Users/richardr/.julia/v0.4/Diversity/src/Collection.jl)

---

<a id="type__generalsimilarity.1" class="lexicon_definition"></a>
#### Diversity.GeneralSimilarity [¶](#type__generalsimilarity.1)
### A general matrix-based Similarity subtype

This subtype of Similarity simply holds a matrix with similarities
between individuals.

#### Members:

- `matrix` A two-dimensional matrix representing similarity between
           individuals. By default this will be the identity matrix,
           but will require the number of species to be instantiated.


*source:*
[Diversity/src/Collection.jl:60](file:///Users/richardr/.julia/v0.4/Diversity/src/Collection.jl)

---

<a id="type__onecommunity.1" class="lexicon_definition"></a>
#### Diversity.Onecommunity [¶](#type__onecommunity.1)
### Partition type allowing only one subcommunity


*source:*
[Diversity/src/Collection.jl:115](file:///Users/richardr/.julia/v0.4/Diversity/src/Collection.jl)

---

<a id="type__subcommunity.1" class="lexicon_definition"></a>
#### Diversity.Subcommunity [¶](#type__subcommunity.1)
### Partition type with multiple subccomunities


*source:*
[Diversity/src/Collection.jl:110](file:///Users/richardr/.julia/v0.4/Diversity/src/Collection.jl)

---

<a id="type__unique.1" class="lexicon_definition"></a>
#### Diversity.Unique [¶](#type__unique.1)
### A subtype of Similarity where all individuals are completely distinct

This type is the simplest Similarity subtype, which identifies all
individuals as unique and completely distinct from each other.


*source:*
[Diversity/src/Collection.jl:16](file:///Users/richardr/.julia/v0.4/Diversity/src/Collection.jl)

## Typealiass [Exported]

---

<a id="typealias__species.1" class="lexicon_definition"></a>
#### Species [¶](#typealias__species.1)
### A subtype of Similarity where all species are completely distinct

This type is the simplest Similarity subtype, which identifies all
species as unique and completely distinct from each other.


*source:*
[Diversity/src/Collection.jl:25](file:///Users/richardr/.julia/v0.4/Diversity/src/Collection.jl)


## Methods [Internal]

---

<a id="method__call.1" class="lexicon_definition"></a>
#### call(::Type{Diversity.GeneralSimilarity},  Z::Array{Float64, 2}) [¶](#method__call.1)
### Constructor for GeneralSimilarity

Creates an instance of the GeneralSimilarity class, with an arbitrary similarity matrix.

#### Arguments:
- `Z`: similarity matrix


*source:*
[Diversity/src/Collection.jl:74](file:///Users/richardr/.julia/v0.4/Diversity/src/Collection.jl)

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
[Diversity/src/CommunityContributions.jl:31](file:///Users/richardr/.julia/v0.4/Diversity/src/CommunityContributions.jl)

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
[Diversity/src/CommunityContributions.jl:31](file:///Users/richardr/.julia/v0.4/Diversity/src/CommunityContributions.jl)

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
[Diversity/src/CommunityContributions.jl:31](file:///Users/richardr/.julia/v0.4/Diversity/src/CommunityContributions.jl)

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
[Diversity/src/CommunityContributions.jl:31](file:///Users/richardr/.julia/v0.4/Diversity/src/CommunityContributions.jl)

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
[Diversity/src/CommunityContributions.jl:31](file:///Users/richardr/.julia/v0.4/Diversity/src/CommunityContributions.jl)

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
[Diversity/src/CommunityContributions.jl:31](file:///Users/richardr/.julia/v0.4/Diversity/src/CommunityContributions.jl)

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
[Diversity/src/EffectiveNumbers.jl:16](file:///Users/richardr/.julia/v0.4/Diversity/src/EffectiveNumbers.jl)

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
[Diversity/src/EffectiveNumbers.jl:16](file:///Users/richardr/.julia/v0.4/Diversity/src/EffectiveNumbers.jl)

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
[Diversity/src/EffectiveNumbers.jl:16](file:///Users/richardr/.julia/v0.4/Diversity/src/EffectiveNumbers.jl)

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
[Diversity/src/Collection.jl:105](file:///Users/richardr/.julia/v0.4/Diversity/src/Collection.jl)

---

<a id="type__similarity.1" class="lexicon_definition"></a>
#### Diversity.Similarity [¶](#type__similarity.1)
### Abstract Similarity supertype for all similarity measures

This type is the abstract superclass of all similarity types. Its
subtypes allow you to define how similarity is measured between
individuals.


*source:*
[Diversity/src/Collection.jl:8](file:///Users/richardr/.julia/v0.4/Diversity/src/Collection.jl)

