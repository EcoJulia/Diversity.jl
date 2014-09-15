## diversity() - calculates subcommunity and ecosystem diversities
##
## Calculates any diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
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
##                     subcommunity diversities
## - returnweights   - boolean describing whether to return subcommunity weights
##
## Returns:
## - some or all (as tuple) of:
##   - vector of ecosystem diversities representing values of q
##   - array of diversities, first dimension representing subcommunities, and
##     last representing values of q
##   - vector of subcommunity weights
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
    
    ## We'll definitely need to calculate subcommunity diversity first
    cd = measure(proportions, qs, Z)

    ## But do we need to calculate anything else?
    if (returnecosystem || returnweights)
        w = mapslices(sum, proportions, 1)
        if (returnecosystem)
            ed = zeros(powers)
            for (i, power) in enumerate(powers)
                ed[i] = powermean(reshape(cd[i, :], size(proportions, 2)),
                                  power, reshape(w, size(proportions, 2)))
            end
            # must be returning ecosystem, but what else?
            return (returncommunity ?
                    (returnweights ? (ed, cd, w) : (ed, cd)) :
                    (returnweights ? (ed, w) : (ed)))
        else # must be returning subcommunity and weights
            return (cd, w)
        end
    else
        # must just be returning subcommunity
        return (cd)
    end
end

## ᾱ() - Normalised similarity-sensitive subcommunity alpha diversity.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
## a vector of qs
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing subcommunities, and
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

subcommunityalphabar = ᾱ

## α() - Raw similarity-sensitive subcommunity alpha diversity.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing subcommunities, and
##   last representing values of q
function α{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("α: Similarity matrix size does not match species number")
    mapslices(p -> powermean(Z * p, qs - 1., p) .^ -1,  proportions, 1)
end

subcommunityalpha = α

## A() - Raw similarity-sensitive ecosystem alpha diversity.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
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
## independent subcommunity counts, for a series of orders, repesented as
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

## ϵ or ρ̄() - Normalised similarity-sensitive subcommunity beta diversity.
##
## β̄ is retained for compatibility (= 1 / ϵ), but we believe ϵ (or ρ̄) to
## be the more fundamental measure.  This is the evenness of the
## subcommunity.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing subcommunities, and
##   last representing values of q
function ϵ{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("β̄: Similarity matrix size does not match species number")

    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions, 1))) / sum(proportions)
    mapslices(p -> powermean(Zp ./ (Z * p), 1. - qs, p),
              proportions * diagm(reshape(mapslices(v -> 1. / sum(v),
                                                    proportions, 1),
                                          (size(proportions, 2)))), 1)
end
subcommunityevenness = ϵ
subcommunityrhobar = ρ̄ = ϵ
function β̄{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    1. ./ ϵ(proportions, qs, Z)
end
subcommunitybetabar = β̄

## ρ() - Raw similarity-sensitive subcommunity beta diversity.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing subcommunities, and
##   last representing values of q
function ρ{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("β: Similarity matrix size does not match species number")

    Zp = Z * reshape(mapslices(sum, proportions, 2), (size(proportions, 1)))
    mapslices(p -> powermean(Zp ./ (Z * p), 1. - qs, p), proportions, 1)
end
subcommunityredundancy = subcommunityrho = ρ
    
function β{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    1. ./ ρ(proportions, qs, Z)
end
subcommunitybeta = β

## R() - Raw similarity-sensitive ecosystem beta diversity.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
R{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(ρ, proportions, qs, Z, true, false, false)
ecosystemredundancy = ecosystemR = R
    
function B{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    1. ./ R(proportions, qs, Z)
end
ecosystemB = B

## E() - Normalised similarity-sensitive ecosystem beta diversity.
##
## B̄ is retained for compatibility (= 1 / E), but we believe E (or R̄) to
## be the more fundamental measure.  This is the average evenness of the
## subcommunities.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - vector of diversities representing values of q
E{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(ϵ, proportions, qs, Z, true, false, false)
ecosystemevenness = E
ecosystemRbar = R̄ = E

function B̄{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    1. ./ R̄(proportions, qs, Z)
end
ecosystemBbar = B̄

## γ̄() - Normalised similarity-sensitive subcommunity gamma diversity.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing subcommunities, and
##   last representing values of q
function γ̄{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("γ̄: Similarity matrix size does not match species number")

    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions, 1))) / sum(proportions)
    mapslices(p -> powermean(Zp, qs - 1., p) .^ -1, proportions, 1)
end

subcommunitygammabar = γ̄

## γ() - Raw similarity-sensitive subcommunity gamma diversity.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
## a vector of qs (or a single number)
##
## Arguments:
## - proportions - population proportions
## - qs - single number or vector of values of parameter q
## - Z - similarity matrix
##
## Returns:
## - array of diversities, first dimension representing subcommunities, and
##   last representing values of q
function γ{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("γ: Similarity matrix size does not match species number")

    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions, 1)))
    mapslices(p -> powermean(Zp, qs - 1., p) .^ -1, proportions, 1)
end

subcommunitygamma = γ

## G() - Raw similarity-sensitive ecosystem gamma diversity.
##
## Calculates diversity of a series of columns representing
## independent subcommunity counts, for a series of orders, repesented as
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
## independent subcommunity counts, for a series of orders, repesented as
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
