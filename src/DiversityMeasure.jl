using Compat

"""
### Enumeration of levels that can exist / be calculated for a metacommunity.
"""
@enum DiversityLevel individualDiversity subcommunityDiversity communityDiversity typeDiversity typeCollectionDiversity metacommunityDiversity

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
abstract DiversityMeasure{FP <: AbstractFloat, AbArray <: AbstractArray}

"""
### Return the ASCII name of the DiversityMeasure

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- String containing simple ASCII name of DiversityMeasure
"""
function getASCIIName(dm::DiversityMeasure)
    s = replace(string(typeof(dm)), "Diversity.", "")
    replace(s, r"{.*}$", "")
end

"""
### Return the character corresponding to the DiversityMeasure

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- String containing unicode (greek) name of DiversityMeasure
"""
function getName
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

@compat (dl::DiversityLevel)(dm::DiversityMeasure) = getPartitionFunction(dm, dl)
@compat (dl::DiversityLevel)(dm::DiversityMeasure, qs) = getPartitionFunction(dm, dl)(qs)

"""
### Metatype of all power mean-based diversity measures

This abstract DiversityMeasure subtype is the metatype of all
diversity measures which are straight power means. PowerMeanMeasure
subtypes allow you to calculate and cache any kind of diversity of a
metacommunity.
"""
abstract PowerMeanMeasure{FP, AbArray} <: DiversityMeasure{FP, AbArray}

"""
### Metatype of all relative entropy-based diversity measures

This abstract DiversityMeasure subtype is the metatype of all
diversity measures which are straight power means.
RelativeEntropyMeasure subtypes allow you to calculate and cache any
kind of diversity of a metacommunity.
"""
abstract RelativeEntropyMeasure{FP, AbArray} <: DiversityMeasure{FP, AbArray}

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

@inline function inddiv(measure::DiversityMeasure, ::Real)
    measure.diversities
end

@inline function inddiv(measure::DiversityMeasure, qs::AbstractVector)
    map(q -> measure.diversities, qs)
end

@inline function inddiv(meta::AbstractMetacommunity, qs)
    map(dm -> inddiv(dm(meta), qs), [RawAlpha, NormalisedAlpha,
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

@inline function subdiv(measure::PowerMeanMeasure, q::Real)
    powermean(inddiv(measure, q), 1.0 - q, measure.abundances)
end

@inline function subdiv(measure::PowerMeanMeasure, qs::AbstractVector)
    map(q -> powermean(inddiv(measure, q),
                       1.0 - q, measure.abundances), qs)
end

@inline function subdiv(measure::RelativeEntropyMeasure, q::Real)
    powermean(inddiv(measure, q), q - 1.0, measure.abundances)
end

@inline function subdiv(measure::RelativeEntropyMeasure, qs::AbstractVector)
    map(q -> powermean(inddiv(measure, q),
                       q - 1.0, measure.abundances), qs)
end

@inline function subdiv(meta::AbstractMetacommunity, qs)
    map(dm -> subdiv(dm(meta), qs), [RawAlpha, NormalisedAlpha,
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

@inline function metadiv(measure::DiversityMeasure, q::Real)
    powermean(vectorise(subdiv(measure, q)), 1.0 - q, measure.weights)
end

@inline function metadiv(measure::DiversityMeasure, qs::AbstractVector)
    map(q -> powermean(vectorise(subdiv(measure, q)),
                       1.0 - q, measure.weights), qs)
end

@inline function metadiv(meta::AbstractMetacommunity, qs)
    map(dm -> metadiv(dm(meta), qs), [RawAlpha, NormalisedAlpha,
                                      RawBeta, NormalisedBeta,
                                      RawRho, NormalisedRho, Gamma])
end

function getPartitionFunction(measure::DiversityMeasure, level::DiversityLevel)
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

- `meta`: a Metacommunity
"""
type RawAlpha{FP, AbArray} <: PowerMeanMeasure{FP, AbArray}
    abundances::AbArray
    weights::Vector{FP}
    diversities::AbArray
end

function RawAlpha(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    w = vectorise(getweight(meta))
    RawAlpha{eltype(ab), typeof(ab)}(ab, w, getordinariness!(meta) .^ -1)
end

getName(::RawAlpha) = "α"
getFullName(::RawAlpha) = "raw alpha diversity"

"""
### Normalised alpha diversity type (ᾱ)

Calculates normalised alpha diversity (ᾱ) of all of the individuals in
a metacommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
type NormalisedAlpha{FP, AbArray} <: PowerMeanMeasure{FP, AbArray}
    abundances::AbArray
    weights::Vector{FP}
    diversities::AbArray
end

function NormalisedAlpha(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    ws = getweight(meta)
    w = vectorise(ws)
    NormalisedAlpha{eltype(ab), typeof(ab)}(ab, w,
                                            ws ./ getordinariness!(meta))
end

getName(::NormalisedAlpha) = "ᾱ"
getFullName(::NormalisedAlpha) = "normalised alpha diversity"

"""
### Distinctiveness (β, raw beta diversity) type

Calculates distinctiveness (β, raw beta diversity) of all of the individuals in a
metacommunity, and caches them for subsequent analysis. This is a
subtype of RelativeEntropyMeasure, meaning that subcommunity and type
composite diversity measures are relative entropies, and their
composite types are powermeans of those measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
type RawBeta{FP, AbArray} <: RelativeEntropyMeasure{FP, AbArray}
    abundances::AbArray
    weights::Vector{FP}
    diversities::AbArray
end

function RawBeta(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    w = vectorise(getweight(meta))
    RawBeta{eltype(ab), typeof(ab)}(ab, w,
                                    getordinariness!(meta) ./
                                    getmetaordinariness!(meta))
end

typealias Distinctiveness RawBeta

getName(::RawBeta) = "β"
getFullName(::RawBeta) = "distinctiveness"

"""
### Normalised beta diversity type (β̄)

Calculates normalised beta diversity (β̄) of all of the individuals in
a metacommunity, and caches them for subsequent analysis. This is a
subtype of RelativeEntropyMeasure, meaning that subcommunity and type
composite diversity measures are relative entropies, and their
composite types are powermeans of those measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
type NormalisedBeta{FP, AbArray} <: RelativeEntropyMeasure{FP, AbArray}
    abundances::AbArray
    weights::Vector{FP}
    diversities::AbArray
end

function NormalisedBeta(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    ws = getweight(meta)
    w = vectorise(ws)
    NormalisedBeta{eltype(ab), typeof(ab)}(ab, w,
                                           getordinariness!(meta) ./
                                           (getmetaordinariness!(meta) .* ws))
end

getName(::NormalisedBeta) = "β̄"
getFullName(::NormalisedBeta) = "effective number of subcommunities"

"""
### Redundancy (ρ, raw beta diversity) type

Calculates redundancy (ρ, raw beta diversity) of all of the
individuals in a metacommunity, and caches them for subsequent
analysis. This is a subtype of PowerMeanMeasure, meaning that all
composite diversity measures are simple powermeans of the individual
measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
type RawRho{FP, AbArray} <: PowerMeanMeasure{FP, AbArray}
    abundances::AbArray
    weights::Vector{FP}
    diversities::AbArray
end

function RawRho(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    w = vectorise(getweight(meta))
    RawRho{eltype(ab), typeof(ab)}(ab, w,
                                   getmetaordinariness!(meta) ./
                                   getordinariness!(meta))
end

typealias Redundancy RawRho

getName(::RawRho) = "ρ"
getFullName(::RawRho) = "redundancy"

"""
### Representativeness (ρ̄, normalised beta diversity) type

Calculates redundancy (ρ̄, normalised beta diversity) of all of the
individuals in a metacommunity, and caches them for subsequent
analysis. This is a subtype of PowerMeanMeasure, meaning that all
composite diversity measures are simple powermeans of the individual
measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
type NormalisedRho{FP, AbArray} <: PowerMeanMeasure{FP, AbArray}
    abundances::AbArray
    weights::Vector{FP}
    diversities::AbArray
end

function NormalisedRho(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    ws = getweight(meta)
    w = vectorise(ws)
    NormalisedRho{eltype(ab), typeof(ab)}(ab, w,
                                          (getmetaordinariness!(meta) .* ws) ./
                                          getordinariness!(meta))
end

typealias Representativeness NormalisedRho

getName(::NormalisedRho) = "ρ̄"
getFullName(::NormalisedRho) = "representativeness"

"""
### Gamma diversity type (γ)

Calculates gamma diversity (γ) of all of the individuals in a
metacommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
type Gamma{FP, AbArray} <: PowerMeanMeasure{FP, AbArray}
    abundances::AbArray
    weights::Vector{FP}
    diversities::AbArray
end

function Gamma(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    ws = getweight(meta)
    w = vectorise(ws)
    Gamma{eltype(ab), typeof(ab)}(ab, w,
                                  ones(eltype(ws), size(ws)) ./
                                  getmetaordinariness!(meta))
end

getName(::Gamma) = "γ"
getFullName(::Gamma) = "gamma diversity"
