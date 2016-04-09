"""
### Calculates subcommunity and supercommunity diversities

Calculates any diversity of a series of columns representing
independent subcommunity counts, for a series of orders, repesented as
a vector of qs, with similarity matrix Z, by default the (naïve)
identity matrix.

#### Arguments:
- `measure`: the diversity function to be used - one of Dα, Dᾱ, Dρ, Dϵ
          (or Dρ̄), Dγ or Dγ̄

- `proportions`:population proportions

- `qs`: single number or vector of values of parameter q

- `Z`: similarity matrix

- `returnsupercommunity`: boolean describing whether to return the
                         supercommunity diversity

- `returnsubcommunity`: boolean describing whether to return the
                       subcommunity diversities 

- `returnweights`: boolean describing whether to return subcommunity weights

#### Returns:
Some or all (as tuple) of:  

- vector of supercommunity diversities representing values of q  

- array of diversities, first dimension representing subcommunities, and
  last representing values of q  

- multidimensional array with dimensions matiching shape of proportions,
  with extra dimension for values of q
"""
function diversity{S <: AbstractFloat}(measure::Function,
                                       proportions::Matrix{S},
                                       qs,
                                       Z::Matrix{S} = eye(size(proportions, 1)),
                                       returnsupercommunity::Bool = true,
                                       returnsubcommunity::Bool = true,
                                       returnweights::Bool = true)
    ## Make sure we actually want to calculate the diversity before
    ## going any further!
    if (!returnsupercommunity && !returnsubcommunity)
        return returnweights ? mapslices(sum, proportions, 1) : nothing
    end

    ## We need our qs to be a vector of floating points
    powers = 1.0 - vec(collect(S, qs))
    
    ## We'll definitely need to calculate subcommunity diversity first
    cd = measure(proportions, qs, Z)

    ## But do we need to calculate anything else?
    if (returnsupercommunity || returnweights)
        w = mapslices(sum, proportions, 1)
        if (returnsupercommunity)
            ed = zeros(powers)
            for (i, power) in enumerate(powers)
                ed[i] = powermean(reshape(cd[i, :], size(proportions, 2)),
                                  power, reshape(w, size(proportions, 2)))
            end
            # must be returning supercommunity, but what else?
            return (returnsubcommunity ?
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
### Raw similarity-sensitive subcommunity alpha diversity / naive-community diversity

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q
"""
subcommunityalpha{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((α, :sub), proportions, qs, sim)

"""
### Normalised similarity-sensitive subcommunity alpha diversity)

Calculates (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q
"""
subcommunityalphabar{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((ᾱ, :sub), proportions, qs, sim)

"""
### Raw similarity-sensitive supercommunity alpha diversity / naive-community diversity

Calculates average raw alpha diversity / naive-community diversity of
a series of subcommunities represented by columns of independent
subcommunity counts, for a series of orders, represented as a vector
of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q
"""
supercommunityA{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((α, :super), proportions, qs, sim)

"""
### Normalised similarity-sensitive supercommunity alpha diversity

Calculates average (normalised alpha) diversity of a series of
subcommunities represented by columns of independent subcommunity
counts, for a series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q
"""
supercommunityAbar{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((ᾱ, :super), proportions, qs, sim)

"""
### Raw similarity-sensitive subcommunity redundancy

Calculates redundancy of a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of redundancies, first dimension representing subcommunities, and
  last representing values of q
"""
subcommunityrho{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((ρ, :sub), proportions, qs, sim)

subcommunityredundancy = subcommunityrho
    
"""
### Raw similarity-sensitive subcommunity beta diversity / distinctiveness / concentration

Calculates the raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q
"""
subcommunitybeta{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    1.0 ./ diversity((ρ, :sub), proportions, -qs, sim)

subcommunitydistinctiveness = subcommunityconcentration = subcommunitybeta

"""
### Normalised similarity-sensitive subcommunity representativeness

Calculates the representativeness of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs. Representativeness
reflects what proportion of the supercommunity each subcommunity is
representative of on average, so if each subcommunity contains 1/xth
of the species, then the average representativeness of the
subcommunities is 1/x.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of representativenesses, first dimension representing subcommunities, and
  last representing values of q
"""
subcommunityrhobar{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((ρ̄, :sub), proportions, qs, sim)

subcommunityrepresentativeness = subcommunityrhobar

"""
### Normalised similarity-sensitive subcommunity beta diversity

Calculates normalised beta diversities or the effective number of
distinct subcommunities perceived by a series of subcommunities
represented by columns of independent subcommunity counts, represented
as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q
"""
subcommunitybetabar{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    1.0 ./ diversity((ρ̄, :sub), proportions, qs, sim)

"""
### Raw similarity-sensitive supercommunity redundancy

Calculates average redundancy of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of redundancies representing values of q
"""
supercommunityR{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((ρ, :super), proportions, qs, sim)

supercommunityredundancy = supercommunityR

"""
### Raw similarity-sensitive supercommunity beta diversity / distinctiveness / concentration

Calculates average raw beta diversity / distinctiveness of or
concentration of species in a series of subcommunities represented by
columns of independent subcommunity counts, for a series of orders,
represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q
"""
function supercommunityB{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique())
    sub, weights = diversity((ρ, Set([:sub, :weights])), proportions, qs, sim)
    powermean(1.0 ./ sub, qs, weights)
end

supercommunitydistinctiveness = supercommunityconcentration = supercommunityB

"""
### Normalised similarity-sensitive supercommunity representativeness

Calculates average representativeness of a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs. Representativeness
reflects what proportion of the supercommunity each subcommunity is
representative of on average, so if each subcommunity contains 1/xth
of the species, then the average representativeness of the
subcommunities is 1/x.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of representativenesses representing values of q
"""
supercommunityRbar{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((ρ̄, :super), proportions, qs, sim)

supercommunityrepresentativeness = supercommunityRbar

"""
### Normalised similarity-sensitive supercommunity beta diversity / effective number of communities

Calculates average normalised beta diversity or the effective number
of distinct subcommunities present in a series of subcommunities
represented by columns of independent subcommunity counts, for a
series of orders, represented as a vector of qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q
"""
function supercommunityBbar{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique())
    sub, weights = diversity((ρ̄, Set([:sub, :weights])), proportions, qs, sim)
    powermean(1.0 ./ sub, qs, weights)
end

"""
### Raw similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q
"""
subcommunitygamma{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((γ, :super), proportions, qs, sim)

"""
### Normalised similarity-sensitive subcommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- array of diversities, first dimension representing subcommunities, and
  last representing values of q
"""
subcommunitygammabar{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((γ̄, :sub), proportions, qs, sim)

"""
### Raw similarity-sensitive supercommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q
"""
supercommunityG{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((γ, :super), proportions, qs, sim)

"""
### Normalised similarity-sensitive supercommunity gamma diversity

Calculates diversity of a series of columns representing independent
subcommunity counts, for a series of orders, represented as a vector of
qs.

#### Arguments:

- `proportions`: population proportions

- `qs`: single number or vector of values of parameter q

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:

- vector of diversities representing values of q
"""
supercommunityGbar{S <: AbstractFloat}(proportions::Matrix{S}, qs, sim = Unique()) =
    diversity((γ̄, :super), proportions, qs, sim)
