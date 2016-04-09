"""
### Calculates subcommunity and supercommunity diversities

Calculates any diversity of a series of columns representing
independent subcommunity counts, for a series of orders, repesented as
a vector of qs, with similarity matrix Z, by default the (naïve)
identity matrix.

#### Arguments:
- Tuple containing:
  * `measure`: the diversity function to be used - one of Diversity.α, Diversity.ᾱ,
Diversity.ρ, Diversity.ρ̄, Diversity.γ or Diversity.γ̄
  * a Set of symbols, containing some or all of :super, :sub and :weights,
detailing what should be returned
- `proportions`:population proportions
- `qs`: single number or vector of values of parameter q
- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

#### Returns:
Some or all (as tuple) of:  

- vector of supercommunity diversities representing values of q  

- array of diversities, first dimension representing subcommunities, and
  last representing values of q  

- multidimensional array with dimensions matiching shape of proportions,
  with extra dimension for values of q
"""
function diversity(kind::@compat(Tuple{Function, Set{Symbol}}),
                   sup::AbstractSupercommunity, qs)
    ## Make sure we actually want to calculate the diversity before
    ## going any further!
    measure, returns = kind
    returnsuper = :supercommunity ∈ returns
    returnsub = :subcommunity ∈ returns
    returnall = :all ∈ returns
    returnweights = :weights ∈ returns
    if (!returnsuper && !returnsub && !returnall)
        return returnweights ? mapslices(sum, getAbundances(sup), 1) : nothing
    end

    ## We need our qs to be a vector of floating points
    powers = 1.0 - vec(collect(S, qs))
    
    ## We'll need to calculate all of the values first
    alld = measure(sup)

    ## Then the subcommunity diversities
    sizes = collect(size(proportions))
    nsub = sizes[2]
    sizes = (length(powers), sizes[2])
    subd = Array(S, sizes)
    for i in 1:nsub
        subd[:,i] = powermean(alld[:,i], powers, proportions[:,i])
    end
    
    ## But do we need to calculate anything else?
    if (returnsuper || returnweights)
        w = mapslices(sum, proportions, 1)
        if (returnsuper)
            superd = zeros(powers)
            for (i, power) in enumerate(powers)
                superd[i] = powermean(reshape(subd[i, :], size(proportions, 2)),
                                  power, reshape(w, size(proportions, 2)))
            end
            # must be returning supercommunity, but what else?
            return (returnsub ?
                    (returnweights ? (superd, subd, w) : (superd, subd)) :
                    (returnweights ? (superd, w) : (superd)))
        else # must be returning subcommunity and weights
            return (subd, w)
        end
    else
        # must just be returning subcommunity
        return (subd)
    end
end

diversity{S <: AbstractFloat}(kind::@compat(Tuple{Function, Set{Symbol}}),
                              proportions::Matrix{S}, qs, z::Matrix{S}) =
                                  diversity(kind, proportions, qs,
                                            GeneralSimilarity(z))

diversity{S <: AbstractFloat}(kind::@compat(Tuple{Function, Symbol}),
                              proportions::Matrix{S}, qs, sim) =
                                  diversity((kind[1], Set{Symbol}([kind[2]])),
                                            proportions, qs, sim)

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

- `z::Matrix{S <: AbstractFloat}` or `sim::T <: Similarity`: similarity matrix or object

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
                                       proportions::Matrix{S}, qs, sim,
                                       returnsupercommunity::Bool = true,
                                       returnsubcommunity::Bool = true,
                                       returnweights::Bool = true)
    returns = Set{Symbol}()
    returnsupercommunity && push!(returns, :super)
    returnsubcommunity && push!(returns, :sub)
    returnweights && push!(returns, :weights)
    diversity((measure, returns), proportions, qs, sim)
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
