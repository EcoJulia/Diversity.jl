# Diversity.Ecology

## Exported

---

<a id="method__generalisedjaccard.1" class="lexicon_definition"></a>
#### generalisedjaccard(proportions::Array{T, 2},  qs) [¶](#method__generalisedjaccard.1)
### Calculate a generalised version of the Jaccard index

Calculates a generalisation of the Jaccard index of a series of
columns representing subcommunity counts. This evaluates to is DG / DA
for a series of orders, repesented as a vector of qs (or a single
number).  It also includes a similarity matrix for the species. This
gives measure of the average distinctiveness of the subcommunities.

#### Arguments:
- `proportions`: population proportions
- `qs`: single number or vector of values of parameter q
- `Z`: similarity matrix

#### Returns:
- Jaccard-related distinctivess measures


*source:*
[Diversity/src/Ecology.jl:143](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

---

<a id="method__generalisedjaccard.2" class="lexicon_definition"></a>
#### generalisedjaccard(proportions::Array{T, 2},  qs,  Z::Array{T, 2}) [¶](#method__generalisedjaccard.2)
### Calculate a generalised version of the Jaccard index

Calculates a generalisation of the Jaccard index of a series of
columns representing subcommunity counts. This evaluates to is DG / DA
for a series of orders, repesented as a vector of qs (or a single
number).  It also includes a similarity matrix for the species. This
gives measure of the average distinctiveness of the subcommunities.

#### Arguments:
- `proportions`: population proportions
- `qs`: single number or vector of values of parameter q
- `Z`: similarity matrix

#### Returns:
- Jaccard-related distinctivess measures


*source:*
[Diversity/src/Ecology.jl:143](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

---

<a id="method__generalisedrichness.1" class="lexicon_definition"></a>
#### generalisedrichness{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2}) [¶](#method__generalisedrichness.1)
### Calculate a generalised version of richness

Calculates (species) richness of a series of columns representing
independent subcommunity counts, which is diversity at q = 0 for any
diversity measure (passed as the second argument). It also includes a
similarity matrix for the species

#### Arguments:
- `measure`: diversity measure to use (one of Dα, Dᾱ, Dρ, Dρ̄, Dγ or Dγ̄)

- `proportions`: population proportions

- `Z`: similarity matrix

### Returns:
- diversity (at ecosystem level) or diversities (of subcommunities)


*source:*
[Diversity/src/Ecology.jl:21](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

---

<a id="method__generalisedrichness.2" class="lexicon_definition"></a>
#### generalisedrichness{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  Z::Array{S<:AbstractFloat, 2}) [¶](#method__generalisedrichness.2)
### Calculate a generalised version of richness

Calculates (species) richness of a series of columns representing
independent subcommunity counts, which is diversity at q = 0 for any
diversity measure (passed as the second argument). It also includes a
similarity matrix for the species

#### Arguments:
- `measure`: diversity measure to use (one of Dα, Dᾱ, Dρ, Dρ̄, Dγ or Dγ̄)

- `proportions`: population proportions

- `Z`: similarity matrix

### Returns:
- diversity (at ecosystem level) or diversities (of subcommunities)


*source:*
[Diversity/src/Ecology.jl:21](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

---

<a id="method__generalisedshannon.1" class="lexicon_definition"></a>
#### generalisedshannon{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2}) [¶](#method__generalisedshannon.1)
### Calculate a generalised version of Shannon entropy

Calculates Shannon entropy of a series of columns representing
independent subcommunity counts, which is log(diversity) at q = 1 for
any diversity measure (passed as the second argument). It also
includes a similarity matrix for the species

#### Arguments:
- `measure`: diversity measure to use (one of Dα, Dᾱ, Dρ, Dρ̄, Dγ or Dγ̄)

- `proportions`: population proportions

- `Z`: similarity matrix

#### Returns:
- entropy (at ecosystem level) or entropies (of subcommunities)


*source:*
[Diversity/src/Ecology.jl:62](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

---

<a id="method__generalisedshannon.2" class="lexicon_definition"></a>
#### generalisedshannon{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  Z::Array{S<:AbstractFloat, 2}) [¶](#method__generalisedshannon.2)
### Calculate a generalised version of Shannon entropy

Calculates Shannon entropy of a series of columns representing
independent subcommunity counts, which is log(diversity) at q = 1 for
any diversity measure (passed as the second argument). It also
includes a similarity matrix for the species

#### Arguments:
- `measure`: diversity measure to use (one of Dα, Dᾱ, Dρ, Dρ̄, Dγ or Dγ̄)

- `proportions`: population proportions

- `Z`: similarity matrix

#### Returns:
- entropy (at ecosystem level) or entropies (of subcommunities)


*source:*
[Diversity/src/Ecology.jl:62](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

---

<a id="method__generalisedsimpson.1" class="lexicon_definition"></a>
#### generalisedsimpson{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2}) [¶](#method__generalisedsimpson.1)
### Calculate a generalised version of Simpson's index

Calculates Simpson's index of a series of columns representing
independent subcommunity counts, which is 1 / diversity at q = 2 for
any diversity measure (passed as the second argument). It also
includes a similarity matrix for the species

#### Arguments:
- `measure`: diversity measure to use (one of Dα, Dᾱ, Dρ, Dρ̄, Dγ or Dγ̄)

- `proportions`: population proportions
- `Z`: similarity matrix

#### Returns:
- concentration (at ecosystem level) or concentrations (of subcommunities)


*source:*
[Diversity/src/Ecology.jl:102](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

---

<a id="method__generalisedsimpson.2" class="lexicon_definition"></a>
#### generalisedsimpson{S<:AbstractFloat}(measure::Function,  proportions::Array{S<:AbstractFloat, 2},  Z::Array{S<:AbstractFloat, 2}) [¶](#method__generalisedsimpson.2)
### Calculate a generalised version of Simpson's index

Calculates Simpson's index of a series of columns representing
independent subcommunity counts, which is 1 / diversity at q = 2 for
any diversity measure (passed as the second argument). It also
includes a similarity matrix for the species

#### Arguments:
- `measure`: diversity measure to use (one of Dα, Dᾱ, Dρ, Dρ̄, Dγ or Dγ̄)

- `proportions`: population proportions
- `Z`: similarity matrix

#### Returns:
- concentration (at ecosystem level) or concentrations (of subcommunities)


*source:*
[Diversity/src/Ecology.jl:102](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

---

<a id="method__jaccard.1" class="lexicon_definition"></a>
#### jaccard{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2}) [¶](#method__jaccard.1)
### Calculate the Jaccard index

Calculates Jaccard index (Jaccard similarity coefficient) of two
columns representing independent subcommunity counts, which is
DA(proportions, 0) / DG(proportions, 0) - 1

#### Arguments:
- `proportions`: population proportions

#### Returns:
- the Jaccard index


*source:*
[Diversity/src/Ecology.jl:163](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

---

<a id="method__richness.1" class="lexicon_definition"></a>
#### richness{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2}) [¶](#method__richness.1)
### Calculate species richness of populations

Calculates (species) richness of a series of columns representing
independent subcommunity counts, which is diversity at q = 0

#### Arguments:
- `proportions`: population proportions

#### Returns:
- diversities of subcommunities


*source:*
[Diversity/src/Ecology.jl:40](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

---

<a id="method__shannon.1" class="lexicon_definition"></a>
#### shannon{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2}) [¶](#method__shannon.1)
### Calculate Shannon entropy of populations

Calculates shannon entropy of a series of columns representing
independent subcommunity counts, which is log(diversity) at q = 1

#### Arguments:
- `proportions`: population proportions

#### Returns:
- entropies of subcommunities


*source:*
[Diversity/src/Ecology.jl:81](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

---

<a id="method__simpson.1" class="lexicon_definition"></a>
#### simpson{S<:AbstractFloat}(proportions::Array{S<:AbstractFloat, 2}) [¶](#method__simpson.1)
### Calculate Simpson's index

Calculates Simpson's index of a series of columns representing
independent subcommunity counts, which is 1 / diversity (or
concentration) at q = 2

#### Arguments:
- `proportions`: population proportions

#### Returns:
- concentrations of subcommunities


*source:*
[Diversity/src/Ecology.jl:122](file:///Users/richardr/.julia/v0.4/Diversity/src/Ecology.jl)

