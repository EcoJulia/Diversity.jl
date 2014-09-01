module Diversity

export powermean, qD, qDZ
export ᾱ, communityalphabar, α, communityalpha, A, ecosystemA, Ā, ecosystemAbar

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
    length(values) == length(weights) ||
    error("Weight and value vectors must be the same length")
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

## qD - calculate Hill number / naive diversity of order q of a
## population with given relative proportions
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - qs - single number or vector of orders of diversity measurement
function qD{S <: FloatingPoint}(proportions::Matrix{S}, qs)
    mapslices((p) ->  qD(p, qs), proportions, 1)
end

## qDZ - calculate Leinster-Cobbold general diversity of order q of a
## population with given relative proportions, and similarity matrix Z
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - q - order of diversity measurement
## - Z - similarity matrix
function qDZ{S <: FloatingPoint,
             T <: Number}(proportions::Vector{S}, q::T,
                          Z::Matrix{S} = eye(length(proportions)))
    l = length(proportions)
    size(Z) == (l, l) ||
    error("Similarity matrix size does not match species number")
    1. / powermean(Z * proportions, q - 1., proportions)
end

## qDZ - calculate general Leinster-Cobbold diversity of >= 1 order q
## of a population with given relative proportions, and similarity
## matrix Z
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - qs - vector of orders of diversity measurement
## - Z - similarity matrix
function qDZ{S <: FloatingPoint,
             T <: Number}(proportions::Vector{S}, qs::Vector{T},
                          Z::Matrix{S} = eye(length(proportions)))
    map((q) ->  qDZ(proportions, q, Z), qs)
end

## qDZ - calculate general Leinster-Cobbold diversity of >= 1 order q
## of a population with given relative proportions, and similarity
## matrix Z
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - qs - single number or vector of orders of diversity measurement
## - Z - similarity matrix
function qDZ{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                                 Z::Matrix{S} = eye(size(proportions)[1]))
    mapslices((p) ->  qDZ(p, qs, Z), proportions, 1)
end

## ᾱ - Normalised similarity-sensitive sub-community alpha diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function ᾱ{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    mapslices((p) ->  qDZ(p, qs, Z),
              proportions * diagm(reshape(mapslices(v -> 1. / sum(v),
                                                    proportions, 1),
                                          (size(proportions)[2]))), 1)
end

communityalphabar = ᾱ

## α - Raw similarity-sensitive sub-community alpha diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function α{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    mapslices((p) ->  qDZ(p, qs, Z),  proportions, 1)
end

communityalpha = α

## A - Raw similarity-sensitive ecosystem alpha diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function A{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    ca = α(proportions, qs, Z)
    indices = 1:length(qs)
    weights = reshape(mapslices(sum, proportions, 1),
                      (size(proportions)[2]))
    map((i) -> powermean(reshape(ca[i, :], (size(proportions)[2])),
                         qs[i], weights), indices)
end

ecosystemA = A

## Ā - Normalised similarity-sensitive ecosystem alpha diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function Ā{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    ca = ᾱ(proportions, qs, Z)
    indices = 1:length(qs)
    weights = reshape(mapslices(sum, proportions, 1), (size(proportions)[2]))
    map((i) -> powermean(reshape(ca[i, :], (size(proportions)[2])),
                         qs[i], weights), indices)
end

ecosystemAbar = Ā

end # module
