"""
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
"""
function powermean{S <: Number,
                   T <: AbstractFloat,
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

## We need to handle lack of automatic promotion between ints and floats in Julia
powermean{T <: Integer,
          U <: Number}(values::Vector{U}, order::T,
                       weights::Vector{U} = ones(values) * 1.) =
                           powermean(values, order * 1., weights)
                           

## We need to handle matrices as well as vectors
function powermean{S <: Number}(values::@compat(Union{Vector{S}, Matrix{S}}),
                                orders::Vector,
                                weights::@compat(Union{Vector, Matrix}) =
                                ones(values) * 1.)
    (size(values) == size(weights)) ||
    error("Values and weights are not the same size")
    map(order -> powermean(values * 1., order * 1., weights * 1.), orders)
end

"""
### Calculates Hill / naive-similarity diversity

Calculates Hill number or naive diversity of order(s) *qs* of a
population with given relative proportions.

#### Arguments:
- `proportions`: relative proportions of different individuals /
               species in population or series of populations
- `qs`: single number or vector of orders of diversity measurement

#### Returns:
- Diversity of order qs (single number or vector of diversities)"""
function qD{S <: AbstractFloat,
            T <: Number}(proportions::Vector{S},
                         qs::@compat(Union{T, Vector{T}}))
    if !isapprox(sum(proportions), 1.)
        warn("qD: Population proportions don't sum to 1, fixing...")
        proportions /= sum(proportions)
    end
    powermean(proportions, qs - 1., proportions) .^ -1
end

## We need to handle matrices as well as vectors
qD{S <: AbstractFloat,
   T <: Number}(proportions::Matrix{S}, qs::@compat(Union{T, Vector{T}})) =
       mapslices(p -> qD(p, qs), proportions, 1)

"""
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
- Diversity of order qs (single number or vector of diversities)"""
function qDZ{S <: AbstractFloat,
             T <: Number}(proportions::Vector{S}, qs::@compat(Union{T, Vector{T}}),
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

## We need to handle matrices as well as vectors
qDZ{S <: AbstractFloat,
    T <: Number}(proportions::Matrix{S}, qs::@compat(Union{T, Vector{T}}),
                 Z::Matrix{S} = eye(size(proportions, 1))) =
                     mapslices(p -> qDZ(p, qs, Z), proportions, 1)
