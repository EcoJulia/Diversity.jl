## diversity() - calculates sub-community and ecosystem diversities
##
## Calculates any diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs, with similarity matrix Z, by default the (naïve) identity
## matrix.
##
## Arguments:
## - measure - the diversity to be used (one of α, ᾱ, β, β̄, γ or γ̄)
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
## - returnecosystem - boolean describing whether to return the
##                     ecosystem diversity
## - returncommunity - boolean describing whether to return the
##                     community diversities
## - returnweights   - boolean describing whether to return community weights
##
## Returns:
## - some or all (as tuple) of:
##   - vector of ecosystem diversities representing values of q
##   - array of diversities, first dimension representing sub-communities, and
##     last representing values of q
##   - vector of community weights
function diversity{S <: FloatingPoint,
                   T <: Number}(measure::Function,
                                proportions::Matrix{S},
                                qs::Union(T, Vector{T}),
                                Z::Matrix{S} = eye(size(proportions, 1)),
                                returnecosystem::Bool = true,
                                returncommunity::Bool = true,
                                returnweights::Bool = true)
    ## Make sure we actually want to calculate the diversity before
    ## going any further!
    if (!returnecosystem && !returncommunity)
        return returnweights ? mapslices(sum, proportions, 1) : nothing
    end

    ## We need our qs to be a vector of floating points
    powers = 1. - convert(Vector{S}, [qs])
    
    ## We'll definitely need to calculate sub-community diversity first
    cd = measure(proportions, qs, Z)

    ## But do we need to calculate anything else?
    if (returnecosystem || returnweights)
        w = mapslices(sum, proportions, 1)
        if (returnecosystem)
            ed = zeros(powers)
            for (i in 1:length(powers))
                ed[i] = powermean(reshape(cd[i, :], size(proportions, 2)),
                                  powers[i], reshape(w, size(proportions, 2)))
            end
            # must be returning ecosystem, but what else?
            return (returncommunity ?
                    (returnweights ? (ed, cd, w) : (ed, cd)) :
                    (returnweights ? (ed, w) : (ed)))
        else # must be returning community and weights
            return (cd, w)
        end
    else
        # must just be returning community
        return (cd)
    end
end

## ᾱ() - Normalised similarity-sensitive sub-community alpha diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
ᾱ{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   mapslices(p -> qDZ(p, qs, Z),
                             proportions *
                             diagm(reshape(mapslices(v -> 1. / sum(v),
                                                     proportions, 1),
                                           (size(proportions, 2)))),
                             1)

communityalphabar = ᾱ

## α() - Raw similarity-sensitive sub-community alpha diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
α{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   mapslices((p) ->  qDZ(p, qs, Z),  proportions, 1)

communityalpha = α

## A() - Raw similarity-sensitive ecosystem alpha diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
A{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(α, proportions, qs, Z, true, false, false)

ecosystemA = A

## Ā() - Normalised similarity-sensitive ecosystem alpha diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
Ā{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(ᾱ, proportions, qs, Z, true, false, false)

ecosystemAbar = Ā

## β̄() - Normalised similarity-sensitive sub-community beta diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function β̄{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions, 1))) / sum(proportions)
    mapslices(p -> powermean((Z * p) ./ Zp, qs - 1., p),
              proportions * diagm(reshape(mapslices(v -> 1. / sum(v),
                                                    proportions, 1),
                                          (size(proportions, 2)))), 1)
end

communitybetabar = β̄

## β() - Raw similarity-sensitive sub-community beta diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function β{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    Zp = Z * reshape(mapslices(sum, proportions, 2), (size(proportions, 1)))
    mapslices(p -> powermean((Z * p) ./ Zp, qs - 1., p), proportions, 1)
end

communitybeta = β

## B() - Raw similarity-sensitive ecosystem beta diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
B{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(β, proportions, qs, Z, true, false, false)

ecosystemB = B

## B̄() - Normalised similarity-sensitive ecosystem beta diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
B̄{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(β̄, proportions, qs, Z, true, false, false)

ecosystemBbar = B̄

## γ̄() - Normalised similarity-sensitive sub-community gamma diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function γ̄{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions, 1))) / sum(proportions)
    mapslices(p -> powermean(Zp, qs - 1., p) .^ -1, proportions, 1)
end

communitygammabar = γ̄

## γ() - Raw similarity-sensitive sub-community gamma diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing sub-communities, and
##   last representing values of q
function γ{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions, 1)))
    mapslices(p -> powermean(Zp, qs - 1., p) .^ -1, proportions, 1)
end

communitygamma = γ

## G() - Raw similarity-sensitive ecosystem gamma diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
G{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(γ, proportions, qs, Z, true, false, false)

ecosystemG = G

## Ḡ() - Normalised similarity-sensitive ecosystem gamma diversity.
##
## Calculates diversity of a series of columns representing
## independent community counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
Ḡ{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(γ̄, proportions, qs, Z, true, false, false)

ecosystemGbar = Ḡ
