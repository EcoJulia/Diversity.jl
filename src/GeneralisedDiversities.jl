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
