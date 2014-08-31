module Diversity

export powermean, qD

## powermean - Calculate order-th power mean of values, weighted by weights
## By default, weights are equal and order is 1, so this is just the arithmetic mean
##
## Arguments:
## - values - values for which to calculate mean
## - order - order of power mean
## - weights - weights of elements, normalised to 1 inside function
##
## Returns:
## - weighted power mean
function powermean{S <: Number, T <: Number, U <: Number}(values::Vector{S},
                   order::T = 1,
                   weights::Vector{U} = ones(FloatingPoint, size(values)))
    ## Normalise weights to sum to 1 (as per Rényi)
    length(values) == length(values) ||
    throw(DimensionMismatch("weight and value vectors must be the same length"))
    proportions = weights / sum(weights)
    power = convert(FloatingPoint, order)
    present = filter(x -> !isapprox(x[1], 0), zip(proportions, values))
    if (isinf(power))
        if (power > 0) # +Inf -> Maximum
            reduce((a, b) -> a[2] > b[2] ? a : b, present)[2]
        else # -Inf -> Minimum
            reduce((a, b) -> a[2] < b[2] ? a : b, present)[2]
        end
    else
        if (isapprox(power, 0))
            mapreduce((pair) -> pair[2] ^ pair[1], *, present)
        else
            mapreduce(pair -> pair[1] * pair[2] ^ power, +,
                      present) ^ (1 / power)
        end
    end
end

## qD - calculate Hill number / naive diversity of order q of a
## population with given relative proportions
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - q - order of diversity measurement
function qD{S <: FloatingPoint, T <: Number}(proportions::Vector{S},
                                             q::T)
  1. / powermean(proportions, q - 1., proportions)
end

## qD - calculate Hill number / naive diversity of order q of a
## population with given relative proportions
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - qs - vector of orders of diversity measurement
function qD{S <: FloatingPoint, T <: Number}(proportions::Vector{S},
                                             qs::Vector{T})
    map((q) ->  qD(proportions, q), qs)
end
end

end # module
