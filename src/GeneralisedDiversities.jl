function ecosystem()
    return [1,2]
end

function subcommunity()
    return [1]
end

function individual()
    return []
end

@doc """
### diversity() - calculates subcommunity and ecosystem diversities

Calculates any diversity of a series of columns representing
independent subcommunity counts, for a series of orders, repesented as
a vector of qs, with similarity matrix Z, by default the (naïve)
identity matrix.

### Arguments:
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

### Returns:
Some or all (as tuple) of:  

* vector of ecosystem diversities representing values of q  

* array of diversities, first dimension representing subcommunities, and
  last representing values of q  

* multidimensional array with dimensions matiching shape of proportions,
  with extra dimension for values of q""" ->
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

@doc """
### Dα() - Raw similarity-sensitive subcommunity alpha diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* array of diversities, first dimension representing subcommunities, and
  last representing values of q""" ->
function Dα{S <: FloatingPoint,
            T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                         Z::Matrix{S} = eye(size(proportions, 1)))
    l = size(proportions, 1)
    size(Z) == (l, l) ||
    error("Dα: Similarity matrix size does not match species number")
    mapslices(p -> powermean(Z * p, qs - 1., p) .^ -1,  proportions, 1)
end

subcommunityalpha = Dα

@doc """
### Dᾱ() - Normalised similarity-sensitive subcommunity alpha diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* array of diversities, first dimension representing subcommunities, and
  last representing values of q""" ->
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

@doc """
### DA() - Raw similarity-sensitive ecosystem alpha diversity

Calculates naive-community diversity of a series of columns
representing independent subcommunity counts, for a series of orders,
represented as a vector of qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* vector of diversities representing values of q""" ->
DA{S <: FloatingPoint,
   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                Z::Matrix{S} = eye(size(proportions, 1))) =
                    diversity(Dα, proportions, qs, Z, true, false, false)
                    
ecosystemA = DA

@doc """
### DĀ() - Normalised similarity-sensitive ecosystem alpha diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* vector of diversities representing values of q""" ->
DĀ{S <: FloatingPoint,
   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                Z::Matrix{S} = eye(size(proportions, 1))) =
                    diversity(Dᾱ, proportions, qs, Z, true, false, false)
                    
ecosystemAbar = DĀ

@doc """
### Dρ() - Raw similarity-sensitive subcommunity beta diversity / redundancy

Calculates redundancy of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* array of redundancies, first dimension representing subcommunities, and
  last representing values of q""" ->
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
    
@doc """
### Dβ() - Raw similarity-sensitive subcommunity beta diversity

Calculates diversities of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

Dβ is retained for compatibility (= 1 / Dρ), but we believe ρ to be the
more fundamental measure - it is the redundancies of the subcommunity.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* array of diversities, first dimension representing subcommunities, and
  last representing values of q""" ->
function β{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    1. ./ ρ(proportions, qs, Z)
end
subcommunitybeta = β

@doc """
### Dϵ or Dρ̄() - Normalised similarity-sensitive subcommunity beta diversity / evenness

Calculates evenness of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* array of evennesses, first dimension representing subcommunities, and
  last representing values of q""" ->
function Dϵ{S <: FloatingPoint,
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
subcommunityevenness = Dϵ
subcommunityrhobar = Dρ̄ = Dϵ


@doc """
### Dβ̄() - Normalised similarity-sensitive subcommunity beta diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

Dβ̄ is retained for compatibility (= 1 / ϵ), but we believe ϵ (or ρ̄) to
be the more fundamental measure - it is the evenness of the
subcommunity.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* array of diversities, first dimension representing subcommunities, and
  last representing values of q""" ->
function Dβ̄{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    1. ./ Dϵ(proportions, qs, Z)
end
subcommunitybetabar = Dβ̄

@doc """
### DR() - Raw similarity-sensitive ecosystem beta diversity / redundancy

Calculates average redundancy of a series of columns representing
independent subcommunity counts, for a series of orders, represented
as a vector of qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* vector of redundancies representing values of q""" ->
DR{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(Dρ, proportions, qs, Z, true, false, false)
ecosystemredundancy = ecosystemR = DR
    
@doc """
### DB() - Raw similarity-sensitive ecosystem beta diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

DB is retained for compatibility (= 1 / DR), but we believe DR to be
the more fundamental measure - it is the average redundancy of the
subcommunities.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* vector of diversities representing values of q""" ->
function DB{S <: FloatingPoint,
           T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                        Z::Matrix{S} = eye(size(proportions, 1)))
    1. ./ DR(proportions, qs, Z)
end
ecosystemB = DB

@doc """
### DE() or DR̄() - Normalised similarity-sensitive ecosystem beta diversity / evenness

Calculates average evenness of a series of columns representing
independent subcommunity counts, for a series of orders, represented
as a vector of qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* vector of evennesses representing values of q""" ->
DE{S <: FloatingPoint,
  T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
               Z::Matrix{S} = eye(size(proportions, 1))) =
                   diversity(Dϵ, proportions, qs, Z, true, false, false)
ecosystemevenness = DE
ecosystemRbar = DR̄ = DE

@doc """
### DB̄() - Normalised similarity-sensitive ecosystem beta diversity

Calculates average diversity of a series of columns representing
independent subcommunity counts, for a series of orders, represented
as a vector of qs.

DB̄ is retained for compatibility, but we believe DE (or DR̄)
to be the more fundamental measure - it is the average evenness of the
subcommunities.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* vector of diversities representing values of q""" ->
function DB̄{S <: FloatingPoint,
            T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                         Z::Matrix{S} = eye(size(proportions, 1)))
    diversity(Dβ̄, proportions, qs, Z, true, false, false)
end
ecosystemBbar = DB̄

@doc """
### Dγ() - Raw similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* array of diversities, first dimension representing subcommunities, and
  last representing values of q""" ->
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

@doc """
### Dγ̄() - Normalised similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* array of diversities, first dimension representing subcommunities, and
last representing values of q""" ->
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

@doc """
### DG() - Raw similarity-sensitive ecosystem gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* vector of diversities representing values of q""" ->
DG{S <: FloatingPoint,
   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                Z::Matrix{S} = eye(size(proportions, 1))) =
                    diversity(Dγ, proportions, qs, Z, true, false, false)
                    
ecosystemG = DG

@doc """
### DḠ() - Normalised similarity-sensitive ecosystem gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

### Arguments:

*proportions* population proportions

*qs* single number or vector of values of parameter q

*Z* similarity matrix

### Returns:

* vector of diversities representing values of q""" ->
DḠ{S <: FloatingPoint,
   T <: Number}(proportions::Matrix{S}, qs::Union(T, Vector{T}),
                Z::Matrix{S} = eye(size(proportions, 1))) =
                    diversity(Dγ̄, proportions, qs, Z, true, false, false)

ecosystemGbar = DḠ
