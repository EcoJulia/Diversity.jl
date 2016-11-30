using Compat

"""
### Enumeration of levels that can exist / be calculated for a metacommunity.
"""
@enum DiversityLevel individualDiversity subcommunityDiversity communityDiversity typeDiversity typeCollectionDiversity metacommunityDiversity metacommunityDiversity

"""
### Generates the function to calculate individual diversities

Generates the function to calculate individual diversities for a
series of orders, represented as a vector of qs.

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- Function which takes a single number or vector of values of
  parameter q, and returns the individual diversities for those
  values.
"""
individualDiversity

"""
### Generates the function to calculate subcommunity diversity

Generates the function to calculate subcommunity diversity for a
series of orders, represented as a vector of qs.

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- Function which takes a single number or vector of values of
  parameter q, and returns the subcommunity diversities for those values.
"""
subcommunityDiversity

"""
### Generates the function to calculate metacommunity diversity

Generates the function to calculate metacommunity diversity for a
series of orders, represented as a vector of qs.

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- Function which takes a single number or vector of values of
  parameter q, and returns the metacommunity diversities for those
  values.
"""
metacommunityDiversity

"""
### DiversityMeasure metatype for all diversity measure types

This type is the abstract metaclass of all diversity measure types.
DiversityMeasure subtypes allow you to calculate and cache any kind of
diversity of a metacommunity.
"""
abstract DiversityMeasure{FP <: AbstractFloat}

"""
### Return the name of the DiversityMeasure

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- String containing true unicode name of DiversityMeasure
"""
function getName{DM <: DiversityMeasure}(dm::DM)
    s = replace(string(typeof(dm)), "Diversity.", "")
    replace(s, r"{.*}$", "")
end

"""
### Return the ASCII name of the DiversityMeasure

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- String containing simple ASCII name of DiversityMeasure
"""
function getASCIIName
end

"""
### Return the full name of the DiversityMeasure

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- String containing full descriptive name of DiversityMeasure
"""
function getFullName
end

@compat (dl::DiversityLevel){DM <: DiversityMeasure}(dm::DM) = getPartitionFunction(dm, dl)
@compat (dl::DiversityLevel){DM <: DiversityMeasure}(dm::DM, qs) = getPartitionFunction(dm, dl)(qs)

"""
### Metatype of all power mean-based diversity measures

This abstract DiversityMeasure subtype is the metatype of all
diversity measures which are straight power means. PowerMeanMeasure
subtypes allow you to calculate and cache any kind of diversity of a
metacommunity.
"""
abstract PowerMeanMeasure{FP} <: DiversityMeasure{FP}

"""
### Metatype of all relative entropy-based diversity measures

This abstract DiversityMeasure subtype is the metatype of all
diversity measures which are straight power means.
RelativeEntropyMeasure subtypes allow you to calculate and cache any
kind of diversity of a metacommunity.
"""
abstract RelativeEntropyMeasure{FP} <: DiversityMeasure{FP}

"""
### Returns individual diversities of a diversity measure

Takes a diversity measure and single order or vector of orders, and
returns the individual diversities for those values.

#### Arguments:

- `dm`: DiversityMeasure
- `q` / `qs`: a single order or a vector of orders

#### Returns:

- Returns individual diversities of `dm` for a single order `q` or a
  vector of order `qs`.
"""
function inddiv
end

@inline function inddiv{DM <: DiversityMeasure}(measure::DM, ::Real)
    measure.diversities
end

@inline function inddiv{DM <: DiversityMeasure,
    Vec <: AbstractVector}(measure::DM, qs::Vec)
    map(q -> measure.diversities, qs)
end

@inline function inddiv{Sup <: AbstractMetacommunity}(sup::Sup, qs)
    map(dm -> inddiv(dm(sup), qs), [RawAlpha, NormalisedAlpha,
                                    RawBeta, NormalisedBeta,
                                    RawRho, NormalisedRho, Gamma])
end

"""
### Calculates subcommunity diversities of a diversity measure

Takes a diversity measure and single order or vector of orders, and
calculates and returns the subcommunity diversities for those values.

#### Arguments:

- `dm`: DiversityMeasure
- `q`: a single order or a vector of orders

#### Returns:

- Returns subcommunity diversities of `dm` for a single order `q` or a
  vector of order `qs`.
"""
function subdiv
end

@inline function subdiv{PM <: PowerMeanMeasure}(measure::PM, q::Real)
    powermean(inddiv(measure, q), 1.0 - q, measure.abundances)
end

@inline function subdiv{PM <: PowerMeanMeasure,
    Vec <: AbstractVector}(measure::PM, qs::Vec)
    map(q -> powermean(inddiv(measure, q),
                       1.0 - q, measure.abundances), qs)
end

@inline function subdiv{RE <: RelativeEntropyMeasure}(measure::RE, q::Real)
    powermean(inddiv(measure, q), q - 1.0, measure.abundances)
end

@inline function subdiv{RE <: RelativeEntropyMeasure,
    Vec <: AbstractVector}(measure::RE, qs::Vec)
    map(q -> powermean(inddiv(measure, q),
                       q - 1.0, measure.abundances), qs)
end

@inline function subdiv{Sup <: AbstractMetacommunity}(sup::Sup, qs)
    map(dm -> subdiv(dm(sup), qs), [RawAlpha, NormalisedAlpha,
                                    RawBeta, NormalisedBeta,
                                    RawRho, NormalisedRho, Gamma])
end

"""
### Calculates metacommunity diversities of a diversity measure

Takes a diversity measure and single order or vector of orders, and
calculates and returns the metacommunity diversities for those values.

#### Arguments:

- `dm`: DiversityMeasure
- `q`: a single order or a vector of orders

#### Returns:

- Returns metacommunity diversities of `dm` for a single order `q` or a
  vector of order `qs`.
"""
function metadiv
end

@inline function metadiv{DM <: DiversityMeasure}(measure::DM, q::Real)
    powermean(arrayise(subdiv(measure, q)), 1.0 - q, measure.weights)
end

@inline function metadiv{DM <: DiversityMeasure,
    Vec <: AbstractVector}(measure::DM, qs::Vec)
    map(q -> powermean(arrayise(subdiv(measure, q)),
                       1.0 - q, measure.weights), qs)
end

@inline function metadiv{Sup <: AbstractMetacommunity}(sup::Sup, qs)
    map(dm -> metadiv(dm(sup), qs), [RawAlpha, NormalisedAlpha,
                                      RawBeta, NormalisedBeta,
                                      RawRho, NormalisedRho, Gamma])
end

function getPartitionFunction{DM <: DiversityMeasure}(measure::DM,
                                                      level::DiversityLevel)
    if (level == individualDiversity)
        return function (qs)
            inddiv(measure, qs)
        end
    elseif (level == subcommunityDiversity)
        function (qs)
            subdiv(measure, qs)
        end
    elseif (level == metacommunityDiversity)
        function (qs)
            metadiv(measure, qs)
        end
    else
        error("Unrecognised diversity level")
    end
end


"""
### Raw alpha diversity type (α)

Calculates raw alpha diversity (α) of all of the individuals in a
metacommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

#### Constructor arguments:

- `sup`: a Metacommunity
"""
type α{FP, AbArray <: AbstractArray, WArray <: AbstractArray} <:
    PowerMeanMeasure{FP}
    abundances::AbArray
    weights::WArray
    diversities::AbArray
end

function α{Sup <: AbstractMetacommunity}(sup::Sup)
    ab = getabundance(sup)
    w = arrayise(slicedim(getweight(sup), 1, 1))
    α{eltype(ab), typeof(ab), typeof(w)}(ab, w,
                                         getordinariness!(sup) .^ -1)
end

typealias RawAlpha α

getASCIIName(::α) = "raw alpha"
getFullName(::α) = "raw alpha diversity"

"""
### Normalised alpha diversity type (ᾱ)

Calculates normalised alpha diversity (ᾱ) of all of the individuals in
a metacommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

#### Constructor arguments:

- `sup`: a Metacommunity
"""
type ᾱ{FP, AbArray <: AbstractArray, WArray <: AbstractArray} <:
    PowerMeanMeasure{FP}
    abundances::AbArray
    weights::WArray
    diversities::AbArray
end

function ᾱ{Sup <: AbstractMetacommunity}(sup::Sup)
    ab = getabundance(sup)
    ws = getweight(sup)
    w = arrayise(slicedim(ws, 1, 1))
    ᾱ{eltype(ab), typeof(ab), typeof(w)}(ab, w,
                                         ws ./ getordinariness!(sup))
end

typealias NormalisedAlpha ᾱ

getASCIIName(::ᾱ) = "normalised alpha"
getFullName(::ᾱ) = "normalised alpha diversity"

"""
### Distinctiveness (β, raw beta diversity) type

Calculates distinctiveness (β, raw beta diversity) of all of the individuals in a
metacommunity, and caches them for subsequent analysis. This is a
subtype of RelativeEntropyMeasure, meaning that subcommunity and type
composite diversity measures are relative entropies, and their
composite types are powermeans of those measures.

#### Constructor arguments:

- `sup`: a Metacommunity
"""
type β{FP, AbArray <: AbstractArray, WArray <: AbstractArray} <:
    RelativeEntropyMeasure{FP}
    abundances::AbArray
    weights::WArray
    diversities::AbArray
end

function β{Sup <: AbstractMetacommunity}(sup::Sup)
    ab = getabundance(sup)
    w = arrayise(slicedim(getweight(sup), 1, 1))
    β{eltype(ab), typeof(ab), typeof(w)}(ab, w,
                                         getordinariness!(sup) ./
                                         getmetaordinariness!(sup))
end

typealias RawBeta β
typealias Distinctiveness β

getASCIIName(::β) = "raw beta"
getFullName(::β) = "distinctiveness"

"""
### Normalised beta diversity type (β̄)

Calculates normalised beta diversity (β̄) of all of the individuals in
a metacommunity, and caches them for subsequent analysis. This is a
subtype of RelativeEntropyMeasure, meaning that subcommunity and type
composite diversity measures are relative entropies, and their
composite types are powermeans of those measures.

#### Constructor arguments:

- `sup`: a Metacommunity
"""
type β̄{FP, AbArray <: AbstractArray, WArray <: AbstractArray} <:
    RelativeEntropyMeasure{FP}
    abundances::AbArray
    weights::WArray
    diversities::AbArray
end

function β̄{Sup <: AbstractMetacommunity}(sup::Sup)
    ab = getabundance(sup)
    ws = getweight(sup)
    w = arrayise(slicedim(ws, 1, 1))
    β̄{eltype(ab), typeof(ab), typeof(w)}(ab, w,
                                         getordinariness!(sup) ./
                                         (getmetaordinariness!(sup) .* ws))
end

typealias NormalisedBeta β̄

getASCIIName(::β̄) = "normalised beta"
getFullName(::β̄) = "effective number of subcommunities"

"""
### Redundancy (ρ, raw beta diversity) type

Calculates redundancy (ρ, raw beta diversity) of all of the
individuals in a metacommunity, and caches them for subsequent
analysis. This is a subtype of PowerMeanMeasure, meaning that all
composite diversity measures are simple powermeans of the individual
measures.

#### Constructor arguments:

- `sup`: a Metacommunity
"""
type ρ{FP, AbArray <: AbstractArray, WArray <: AbstractArray} <:
    PowerMeanMeasure{FP}
    abundances::AbArray
    weights::WArray
    diversities::AbArray
end

function ρ{Sup <: AbstractMetacommunity}(sup::Sup)
    ab = getabundance(sup)
    w = arrayise(slicedim(getweight(sup), 1, 1))
    ρ{eltype(ab), typeof(ab), typeof(w)}(ab, w,
                                         getmetaordinariness!(sup) ./
                                         getordinariness!(sup))
end

typealias RawRho ρ
typealias Redundancy ρ

getASCIIName(::ρ) = "raw rho"
getFullName(::ρ) = "redundancy"

"""
### Representativeness (ρ̄, normalised beta diversity) type

Calculates redundancy (ρ̄, normalised beta diversity) of all of the
individuals in a metacommunity, and caches them for subsequent
analysis. This is a subtype of PowerMeanMeasure, meaning that all
composite diversity measures are simple powermeans of the individual
measures.

#### Constructor arguments:

- `sup`: a Metacommunity
"""
type ρ̄{FP, AbArray <: AbstractArray, WArray <: AbstractArray} <:
    PowerMeanMeasure{FP}
    abundances::AbArray
    weights::WArray
    diversities::AbArray
end

function ρ̄{Sup <: AbstractMetacommunity}(sup::Sup)
    ab = getabundance(sup)
    ws = getweight(sup)
    w = arrayise(slicedim(ws, 1, 1))
    ρ̄{eltype(ab), typeof(ab), typeof(w)}(ab, w,
                                         (getmetaordinariness!(sup) .* ws) ./
                                         getordinariness!(sup))
end

typealias NormalisedRho ρ̄
typealias Representativeness ρ̄

getASCIIName(::ρ̄) = "normalised rho"
getFullName(::ρ̄) = "representativeness"

"""
### Gamma diversity type (γ)

Calculates gamma diversity (γ) of all of the individuals in a
metacommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

#### Constructor arguments:

- `sup`: a Metacommunity
"""
type γ{FP, AbArray <: AbstractArray, WArray <: AbstractArray} <:
    PowerMeanMeasure{FP}
    abundances::AbArray
    weights::WArray
    diversities::AbArray
end

function γ{Sup <: AbstractMetacommunity}(sup::Sup)
    ab = getabundance(sup)
    ws = getweight(sup)
    w = arrayise(slicedim(ws, 1, 1))
    γ{eltype(ab), typeof(ab), typeof(w)}(ab, w,
                                         ones(eltype(ws), size(ws)) ./
                                         getmetaordinariness!(sup))
end

typealias Gamma γ

getASCIIName(::γ) = "gamma"
getFullName(::γ) = "gamma diversity"
