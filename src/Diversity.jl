module Diversity

export powermean, qD, qDZ
export ᾱ, communityalphabar, α, communityalpha, A, ecosystemA, Ā, ecosystemAbar
export β̄, communitybetabar, β, communitybeta, B, ecosystemB, B̄, ecosystemBbar
export γ̄, communitygammabar, γ, communitygamma, G, ecosystemG, Ḡ, ecosystemGbar

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

function powermean{S <: Number, T <: Number, U <: Number}(values::Vector{S},
                   orders::Vector{T},
                   weights::Vector{U} = ones(FloatingPoint, size(values)))
    map((order) ->  powermean(values, order, weights), orders)
end

## qD - calculate Hill number / naive diversity of order q of a
## population with given relative proportions
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - qs - single number or vector of orders of diversity measurement
function qD{S <: FloatingPoint}(proportions::Vector{S}, qs)
  powermean(proportions, qs - 1., proportions) .^ -1
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

## qDZ - calculate Leinster-Cobbold general diversity of >= 1 order q
## of a population with given relative proportions, and similarity
## matrix Z
##
## Arguments:
## - proportions - relative proportions of different individuals /
##                 species in population
## - qs - single number or vector of orders of diversity measurement
## - Z - similarity matrix
function qDZ{S <: FloatingPoint}(proportions::Vector{S}, qs,
                                 Z::Matrix{S} = eye(length(proportions)))
    l = length(proportions)
    size(Z) == (l, l) ||
    error("Similarity matrix size does not match species number")
    powermean(Z * proportions, qs - 1., proportions) .^ -1
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
## - proportions - population proportions
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
## - proportions - population proportions
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
## - proportions - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
function A{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    ca = α(proportions, qs, Z)
    indices = 1:length(qs)
    weights = reshape(mapslices(sum, proportions, 1),
                      (size(proportions)[2]))
    map((i) -> powermean(reshape(ca[i, :], (size(proportions)[2])),
                         1 - qs[i], weights), indices)
end

ecosystemA = A

## Ā - Normalised similarity-sensitive ecosystem alpha diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - proportions - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
function Ā{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    ca = ᾱ(proportions, qs, Z)
    indices = 1:length(qs)
    weights = reshape(mapslices(sum, proportions, 1), (size(proportions)[2]))
    map((i) -> powermean(reshape(ca[i, :], (size(proportions)[2])),
                         1 - qs[i], weights), indices)
end

ecosystemAbar = Ā

## β̄ - Normalised similarity-sensitive sub-community beta diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - proportions - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function β̄{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions)[1])) / sum(proportions)
    mapslices((p) ->  powermean((Z * p) ./ Zp, qs - 1., p),
              proportions * diagm(reshape(mapslices(v -> 1. / sum(v),
                                                    proportions, 1),
                                          (size(proportions)[2]))), 1)
end

communitybetabar = β̄

## β - Raw similarity-sensitive sub-community beta diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - proportions - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function β{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    Zp = Z * reshape(mapslices(sum, proportions, 2), (size(proportions)[1]))
    mapslices((p) ->  powermean((Z * p) ./ Zp, qs - 1., p), proportions, 1)
end

communitybeta = β

## B - Raw similarity-sensitive ecosystem beta diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - proportions - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
function B{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    ca = β(proportions, qs, Z)
    indices = 1:length(qs)
    weights = reshape(mapslices(sum, proportions, 1), (size(proportions)[2]))
    map((i) -> powermean(reshape(ca[i, :], (size(proportions)[2])),
                         1 - qs[i], weights), indices)
end

ecosystemB = B

## B̄ - Normalised similarity-sensitive ecosystem beta diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - proportions - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
function B̄{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    ca = β̄(proportions, qs, Z)
    indices = 1:length(qs)
    weights = reshape(mapslices(sum, proportions, 1), (size(proportions)[2]))
    map((i) -> powermean(reshape(ca[i, :], (size(proportions)[2])),
                         1 - qs[i], weights), indices)
end

ecosystemBbar = B̄

## γ̄ - Normalised similarity-sensitive sub-community gamma diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - proportions - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function γ̄{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions)[1])) / sum(proportions)
    mapslices((p) ->  powermean(Zp, qs - 1., p) .^ -1, proportions, 1)
end

communitygammabar = γ̄

## γ - Raw similarity-sensitive sub-community gamma diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - proportions - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function γ{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions)[1]))
    mapslices((p) ->  powermean(Zp, qs - 1., p) .^ -1, proportions, 1)
end

communitygamma = γ

## G - Raw similarity-sensitive ecosystem gamma diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - proportions - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
function G{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    ca = γ(proportions, qs, Z)
    indices = 1:length(qs)
    weights = reshape(mapslices(sum, proportions, 1), (size(proportions)[2]))
    map((i) -> powermean(reshape(ca[i, :], (size(proportions)[2])),
                         1 - qs[i], weights), indices)
end

ecosystemG = G

## Ḡ - Normalised similarity-sensitive ecosystem gamma diversity.
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - proportions - population proportions
## - qs - vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
function Ḡ{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                               Z::Matrix{S} = eye(size(proportions)[1]))
    ca = γ̄(proportions, qs, Z)
    indices = 1:length(qs)
    weights = reshape(mapslices(sum, proportions, 1), (size(proportions)[2]))
    map((i) -> powermean(reshape(ca[i, :], (size(proportions)[2])),
                         1 - qs[i], weights), indices)
end

ecosystemGbar = Ḡ

end # module
