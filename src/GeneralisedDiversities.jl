"""
!!summary(Calculates subcommunity and ecosystem diversities)

Calculates any diversity of a series of columns representing
independent subcommunity counts, for a series of orders, repesented as
a vector of qs, with similarity matrix Z, by default the (naïve)
identity matrix.

#### Arguments:
*measure* the diversity function to be used - one of Dα, Dᾱ, Dρ, Dϵ
          (or Dρ̄), Dγ or Dγ̄

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

*returnecosystem* boolean describing whether to return the
                  ecosystem diversity

*returncommunity* boolean describing whether to return the
                  subcommunity diversities 

*returnweights* boolean describing whether to return subcommunity weights

#### Returns:
Some or all (as tuple) of:  

* vector of ecosystem diversities representing values of q  

* array of diversities, first dimension representing subcommunities, and
  last representing values of q  

* multidimensional array with dimensions matiching shape of proportions,
  with extra dimension for values of q
"""
:diversity

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
    powers = 1. - convert(Vector{S}, collect(qs))
    
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

"""
!!summary(Raw similarity-sensitive subcommunity alpha diversity / naive-community diversity)

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* array of diversities, first dimension representing subcommunities, and
  last representing values of q
"""
[:Dα, :subcommunityalpha]

function Dα{S <: FloatingPoint,
            T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                         Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("Dα: Similarity matrix size does not match species number")
    mapslices(p -> powermean(Z * p, qs - 1., p) .^ -1,  proportions, 1)
end

subcommunityalpha = Dα

"""
!!summary(Normalised similarity-sensitive subcommunity alpha diversity)

Calculates (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* array of diversities, first dimension representing subcommunities, and
  last representing values of q
"""
[:Dᾱ, :subcommunityalphabar]

Dᾱ{S <: FloatingPoint,
   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                Z::Matrix{S} = eye(size(proportions, 1))) =
                    mapslices(p -> qDZ(p, qs, Z),
                              proportions *
                              diagm(reshape(mapslices(v -> 1. / sum(v),
                                                      proportions, 1),
                                            (size(proportions, 2)))),
                              1)
                              
subcommunityalphabar = Dᾱ

"""
!!summary(Raw similarity-sensitive ecosystem alpha diversity / naive-community diversity)

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector
of qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* vector of diversities representing values of q
"""
[:DA, :ecosystemA]

DA{S <: FloatingPoint,
   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                Z::Matrix{S} = eye(size(proportions, 1))) =
                    diversity(Dα, proportions, qs, Z, true, false, false)
                    
ecosystemA = DA

"""
!!summary(Normalised similarity-sensitive ecosystem alpha diversity)

Calculates average (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* vector of diversities representing values of q
"""
[:DĀ, :ecosystemAbar]

DĀ{S <: FloatingPoint,
   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                Z::Matrix{S} = eye(size(proportions, 1))) =
                    diversity(Dᾱ, proportions, qs, Z, true, false, false)
                    
ecosystemAbar = DĀ

"""
!!summary(Raw similarity-sensitive subcommunity redundancy)

Calculates redundancy of a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* array of redundancies, first dimension representing subcommunities, and
  last representing values of q
"""
[:Dρ, :subcommunityrho, :subcommunityredundancy]

function Dρ{S <: FloatingPoint,
            T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                         Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("Dρ: Similarity matrix size does not match species number")
    
    Zp = Z * reshape(mapslices(sum, proportions, 2), (size(proportions, 1)))
    mapslices(p -> powermean(Zp ./ (Z * p), 1. - qs, p), proportions, 1)
end

subcommunityredundancy = subcommunityrho = Dρ
    
"""
!!summary(Raw similarity-sensitive subcommunity beta diversity / distinctiveness / concentration)

Calculates the raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* array of diversities, first dimension representing subcommunities, and
  last representing values of q
"""
[:Dβ, :subcommunitybeta, :subcommunitydistinctiveness,
  :subcommunityconcentration]
  
function Dβ{S <: FloatingPoint,
            T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                         Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("Dβ: Similarity matrix size does not match species number")
    
    Zp = Z * reshape(mapslices(sum, proportions, 2), (size(proportions, 1)))
    mapslices(p -> powermean((Z * p) ./ Zp, qs - 1., p), proportions, 1)
end

subcommunitybeta = subcommunitydistinctiveness =
    subcommunityconcentration = Dβ

"""
!!summary(Normalised similarity-sensitive subcommunity representativeness)

Calculates the representativeness of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs. Representativeness
reflects what proportion of the ecosystem each subcommunity is
representative of on average, so if each subcommunity contains 1/xth
of the species, then the average representativeness of the
subcommunities is 1/x.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* array of representativenesses, first dimension representing subcommunities, and
  last representing values of q
"""
[:Dρ̄, :subcommunityrhobar, :subcommunityrepresentativeness, :Dϵ,
  :subcommunityepsilon]

function Dρ̄{S <: FloatingPoint,
            T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                         Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("Dϵ: Similarity matrix size does not match species number")

    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions, 1))) / sum(proportions)
    mapslices(p -> powermean(Zp ./ (Z * p), 1. - qs, p),
              proportions * diagm(reshape(mapslices(v -> 1. / sum(v),
                                                    proportions, 1),
                                          (size(proportions, 2)))), 1)
end
subcommunityrhobar = subcommunityrepresentativeness = Dρ̄
subcommunityepsilon = Dϵ = Dρ̄

"""
!!summary(Normalised similarity-sensitive subcommunity beta diversity)

Calculates normalised beta diversities or the effective number of
distinct subcommunities perceived by a series of subcommunities
represented by columns of independent subcommunity counts, represented
as a vector of qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* array of diversities, first dimension representing subcommunities, and
  last representing values of q
"""
[:Dβ̄, :subcommunitybetabar]

function Dβ̄{S <: FloatingPoint,
            T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                         Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("Dβ̄: Similarity matrix size does not match species number")
    
    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions, 1))) / sum(proportions)
    mapslices(p -> powermean((Z * p) ./ Zp, qs - 1., p),
              proportions * diagm(reshape(mapslices(v -> 1. / sum(v),
                                                    proportions, 1),
                                          (size(proportions, 2)))), 1)
end
subcommunitybetabar = Dβ̄

"""
!!summary(Raw similarity-sensitive ecosystem redundancy)

Calculates average redundancy of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* vector of redundancies representing values of q
"""
[:DR, :ecosystemR, :ecosystemredundancy]

DR{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(Dρ, proportions, qs, Z, true, false, false)
ecosystemredundancy = ecosystemR = DR
    
"""
!!summary(Raw similarity-sensitive ecosystem beta diversity / distinctiveness / concentration)

Calculates average raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* vector of diversities representing values of q
"""
[:DB, :ecosystemB, :ecosystemdistinctiveness, :ecosystemconcentration]

function DB{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    diversity(Dβ, proportions, qs, Z, true, false, false)
end
ecosystemB = ecosystemdistinctiveness = ecosystemconcentration = DB

"""
!!summary(Normalised similarity-sensitive ecosystem representativeness)

Calculates average representativeness of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs. Representativeness
reflects what proportion of the ecosystem each subcommunity is
representative of on average, so if each subcommunity contains 1/xth
of the species, then the average representativeness of the
subcommunities is 1/x.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* vector of representativenesses representing values of q
"""
[:DR̄, :ecosystemRbar, :ecosystemrepresentativeness, :DE, :ecosystemE]

DR̄{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(Dρ̄, proportions, qs, Z, true, false, false)

ecosystemRbar = ecosystemrepresentativeness = DR̄
ecosystemE = DE = DR̄

"""
!!summary(Normalised similarity-sensitive ecosystem beta diversity / effective number of communities)

Calculates average normalised beta diversity or the effective number
of distinct subcommunities present in a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* vector of diversities representing values of q
"""
[:DB̄, :ecosystemBbar]

function DB̄{S <: FloatingPoint,
            T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                         Z::Matrix{S} = eye(size(proportions, 1)))
    diversity(Dβ̄, proportions, qs, Z, true, false, false)
end
ecosystemBbar = DB̄

"""
!!summary(Raw similarity-sensitive subcommunity gamma diversity)

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* array of diversities, first dimension representing subcommunities, and
  last representing values of q
"""
[:Dγ, :subcommunitygamma]

function Dγ{S <: FloatingPoint,
            T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                         Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("Dγ: Similarity matrix size does not match species number")

    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions, 1)))
    mapslices(p -> powermean(Zp, qs - 1., p) .^ -1, proportions, 1)
end

subcommunitygamma = Dγ

"""
!!summary(Normalised similarity-sensitive subcommunity gamma diversity)

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* array of diversities, first dimension representing subcommunities, and
last representing values of q
"""
[:Dγ̄, :subcommunitygammabar]

function Dγ̄{S <: FloatingPoint,
            T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                         Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("Dγ̄: Similarity matrix size does not match species number")

    Zp = Z * reshape(mapslices(sum, proportions, 2),
                     (size(proportions, 1))) / sum(proportions)
    mapslices(p -> powermean(Zp, qs - 1., p) .^ -1, proportions, 1)
end

subcommunitygammabar = Dγ̄

"""
!summary(Raw similarity-sensitive ecosystem gamma diversity)

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* vector of diversities representing values of q
"""
[:DG, :ecosystemG]

DG{S <: FloatingPoint,
   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                Z::Matrix{S} = eye(size(proportions, 1))) =
                    diversity(Dγ, proportions, qs, Z, true, false, false)

ecosystemG = DG

"""
!!summary(Normalised similarity-sensitive ecosystem gamma diversity)

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

#### Returns:

* vector of diversities representing values of q
"""
[:DḠ, :ecosystemGbar]

DḠ{S <: FloatingPoint,
   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                Z::Matrix{S} = eye(size(proportions, 1))) =
                    diversity(Dγ̄, proportions, qs, Z, true, false, false)

ecosystemGbar = DḠ
