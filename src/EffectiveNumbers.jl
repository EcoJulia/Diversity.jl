Docile.@doc """
### powermean()

Calculates *order*th power mean of *values*, weighted by
*weights*. By default, *weights* are equal and *order*
is 1, so this is just the arithmetic mean.

### Arguments:
* values: values for which to calculate mean
* order: order of power mean
* weights: weights of elements, normalised to 1 inside function

### Returns:
* weighted power mean""" ->
function powermean{S <: Number,
                   T <: FloatingPoint,
                   U <: Number}(values::Vector{S},
                                order::T = 1.,
                                weights::Vector{U} = ones(values) * 1.)
    ## Normalise weights to sum to 1 (as per Rényi)
    length(values) == length(weights) ||
    error("powermean: Weight and value vectors must be the same length")
    proportions = weights / sum(weights)
    present = filter(x -> !isapprox(x[1], 0.), zip(proportions, values))
    if (isinf(order))
        if (order > 0.) # +Inf -> Maximum
            reduce((a, b) -> a[2] > b[2] ? a : b, present)[2]
        else # -Inf -> Minimum
            reduce((a, b) -> a[2] < b[2] ? a : b, present)[2]
        end
    else
        if (isapprox(order, 0))
            mapreduce(pair -> pair[2] ^ pair[1], *, present)
        else
            mapreduce(pair -> pair[1] * pair[2] ^ order, +,
                      present) ^ (1. / order)
        end
    end
end

Docile.@doc """
Calculates *order*th power mean of *values*, weighted by *weights*.
By default *weights* are equal, and *order* is 1, so this is just
the arithmetic mean.""" -> powermean

## We need to handle lack of automatic promotion between ints and floats in Julia
powermean{T <: Integer,
          U <: Number}(values::Vector{U}, order::T,
                       weights::Vector{U} = ones(values) * 1.) =
                           powermean(values, order * 1., weights)
                           
Docile.@doc """
### powermean()

Calculates power mean of *values* of order *orders*, weighted by
*weights*. By default, *weights* are equal and *order* is 1, so this
is just the arithmetic mean.

### Arguments:
* values: values for which to calculate mean
* orders: vector of orders of power mean
* weights: weights of elements, normalised to 1 inside function

### Returns:
* weighted power means""" ->
function powermean{S <: Number}(values::Union(Vector{S}, Matrix{S}),
                                orders::Vector,
                                weights::Union(Vector, Matrix) =
                                ones(values) * 1.)
    (size(values) == size(weights)) ||
    error("Values and weights are not the same size")
    map(order -> powermean(values * 1., order * 1., weights * 1.), orders)
end

Docile.@doc """
### qD()

Calculates Hill number or naive diversity of order(s) *qs* of a
population with given relative proportions.

### Arguments:
* proportions: relative proportions of different individuals /
               species in population
* qs: single number or vector of orders of diversity measurement

### Returns:
* Diversity of order qs (single number or vector of diversities)""" ->
function qD{S <: FloatingPoint,
            T <: Number}(proportions::Vector{S},
                         qs::Union(T, Vector{T}))
    if !isapprox(sum(proportions), 1.)
        warn("qD: Population proportions don't sum to 1, fixing...")
        proportions /= sum(proportions)
    end
    powermean(proportions, qs - 1., proportions) .^ -1
end

Docile.@doc """
### qD()

Calculates Hill number or naive diversity of orders *qs* of a
population or series of subpopulations with given relative
*proportions*.""" -> qD

Docile.@doc """
### qD()

Calculates Hill number or naive diversity of order(s) *qs* of a series
of subpopulations with given relative *proportions*.

### Arguments:
* proportions - relative proportions of different individuals /
                species in a series of subpopulations
* qs - single number or vector of orders of diversity measurement

### Returns:
* Diversities of order qs (single number or vector of diversities)""" ->
qD{S <: FloatingPoint,
   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T})) =
       mapslices(p -> qD(p, qs), proportions, 1)

Docile.@doc """
### qDZ()

Calculates Leinster-Cobbold general diversity of >= 1 order(s) *qs* of
a population with given relative *proportions*, and similarity matrix
*Z*.

### Arguments:
* proportions: relative proportions of different individuals /
               species in population
* qs: single number or vector of orders of diversity measurement
* Z: similarity matrix

### Returns:
* Diversity of order qs (single number or vector of diversities)""" ->
function qDZ{S <: FloatingPoint,
             T <: Number}(proportions::Vector{S}, qs::Union(T, Vector{T}),
                          Z::Matrix{S} = eye(length(proportions)))
    if !isapprox(sum(proportions), 1.)
        warn("qDZ: Population proportions don't sum to 1, fixing...")
        proportions /= sum(proportions)
    end

    l = length(proportions)
    size(Z) == (l, l) ||
    error("qDZ: Similarity matrix size does not match species number")
    powermean(Z * proportions, qs - 1., proportions) .^ -1
end

Docile.@doc """
### qDZ()

Calculates Leinster-Cobbold general diversity of >= 1 order(s) *qs* of
a population of series of subpopulations with given relative
*proportions*, and similarity matrix *Z*.""" -> qDZ

Docile.@doc """
### qDZ()

Calculates Leinster-Cobbold general diversity of >= 1 order *qs* of a
series of subpopulations with given relative *proportions*, and
similarity matrix *Z*.

### Arguments:
* proportions - relative proportions of different individuals /
                species in a series of subpopulations
* qs: single number or vector of orders of diversity measurement
* Z: similarity matrix

### Returns:
* Diversities of order qs (single number or vector of diversities)""" ->
qDZ{S <: FloatingPoint,
    T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                 Z::Matrix{S} = eye(size(proportions, 1))) =
                     mapslices(p -> qDZ(p, qs, Z), proportions, 1)
