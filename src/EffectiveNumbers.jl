## powermean - Calculates order-th power mean of values, weighted by
## weights By default, weights are equal and order is 1, so this is
## just the arithmetic mean
##
## Arguments:
## - values - values for which to calculate mean
## - order - order of power mean
## - weights - weights of elements, normalised to 1 inside function
##
## Returns:
## - weighted power mean
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

## We need to handle lack of automatic promotion between ints and floats in Julia
powermean{T <: Integer,
          U <: Number}(values, order::T,
                       weights::Vector{U} = ones(values) * 1.) =
                           powermean(values, convert(U, order), weights)
                           
## Handle the likelihood of multiple orders of the mean being needed
powermean{S <: Number}(values::Vector{S}, orders::Vector,
                       weights::Vector = ones(values) * 1.) =
                           map(order ->
                               powermean(values, order * 1., weights * 1.),
                               orders)

## qD - calculates Hill number / naive diversity of order q of a
## population with given relative proportions
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - qs - single number or vector of orders of diversity measurement
##
## Returns:
## - Diversity of order qs (single number or vector of diversities)
qD{S <: FloatingPoint,
   T <: Number}(proportions::Vector{S},
                qs::Union(T, Vector{T})) =
                    powermean(proportions, qs - 1., proportions) .^ -1

## qD - calculates Hill number / naive diversity of order q of a
## population with given relative proportions
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - qs - single number or vector of orders of diversity measurement
##
## Returns:
## - Diversity of order qs (single number or vector of diversities)
qD{S <: FloatingPoint,
   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T})) =
       mapslices(p -> qD(p, qs), proportions, 1)
       
## qDZ - calculates Leinster-Cobbold general diversity of >= 1 order q
## of a population with given relative proportions, and similarity
## matrix Z
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - qs - single number or vector of orders of diversity measurement
## - Z - similarity matrix
##
## Returns:
## - Diversity of order qs (single number or vector of diversities)
function qDZ{S <: FloatingPoint,
             T <: Number}(proportions::Vector{S}, qs::Union(T, Vector{T}),
                          Z::Matrix{S} = eye(length(proportions)))
    l = length(proportions)
    size(Z) == (l, l) ||
    error("qDZ: Similarity matrix size does not match species number")
    powermean(Z * proportions, qs - 1., proportions) .^ -1
end

## qDZ - calculates general Leinster-Cobbold diversity of >= 1 order q
## of a population with given relative proportions, and similarity
## matrix Z
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - qs - single number or vector of orders of diversity measurement
## - Z - similarity matrix
##
## Returns:
## - Diversity of order qs (single number or vector of diversities)
qDZ{S <: FloatingPoint,
    T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                 Z::Matrix{S} = eye(size(proportions, 1))) =
                     mapslices(p -> qDZ(p, qs, Z), proportions, 1)
