using Diversity

@doc """
### jostD()

Calculates Hill number of order q of population(s) with given relative
proportions.

### Arguments:
*proportions* relative proportions of different individuals / species
              in population (vector, or matrix where columns are
              individual populations)

*qs* single number or vector of orders of diversity measurement

### Returns:
* Diversity of order qs (single number or vector of diversities)""" ->
function jostD(proportions, qs)
    qD(proportions, qs)
end

@doc """
### jostalpha() (or jostα())

Calculates Jost's alpha diversity of a series of columns representing
independent community counts, for a series of orders, repesented as a
vector of qs. This is just the naive-community ecosystem diversity
divided by the naive-community beta diversity.

### Arguments:
*proportions* relative proportions of different individuals / species
              in population (vector, or matrix where columns are
              for individual sub-communities)

*qs* single number or vector of orders of diversity measurement

### Returns:
* array of diversities, first dimension representing sub-communities, and
  last representing values of q""" ->
function jostalpha{S <: FloatingPoint,
                   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}))
    DA(proportions, qs) ./
    qD(reshape(mapslices(sum, proportions, (1,)), size(proportions)[2]), qs)
end

jostα = jostalpha

@doc """
### jostbeta() (or jostβ())

Calculates Jost's beta diversity of a series of columns representing
independent community counts, for a series of orders, repesented as a
vector of qs. This is just the naive gamma diversity divided by
Jost's alpha diversity

### Arguments:
*proportions* relative proportions of different individuals / species
              in population (vector, or matrix where columns are
              for individual sub-communities)

*qs* single number or vector of orders of diversity measurement

### Returns:
* array of diversities, first dimension representing sub-communities, and
  last representing values of q""" ->
function jostbeta{S <: FloatingPoint,
                  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}))
    DḠ(proportions, qs) ./ jostalpha(proportions, qs)
end

jostβ = jostbeta
