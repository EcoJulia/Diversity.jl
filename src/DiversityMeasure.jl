using Compat
using DataFrames

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
    DiversityMeasure

This type is the abstract supertype of all diversity measure types.
DiversityMeasure subtypes allow you to calculate and cache any kind of
diversity of a metacommunity.
"""
@compat abstract type DiversityMeasure{FP <: AbstractFloat, AbArray <: AbstractArray} end

"""
    getscnames(dm::DiversityMeasure)

Return the names of the subcommunities of the metacommunity being analysed

# Arguments:

- `dm`: DiversityMeasure

# Returns:

- Vector of Strings of names of subcommunities.
"""
function getscnames(dm::DiversityMeasure)
    return dm.scnames
end

"""
    gettypenames(dm::DiversityMeasure)

Return the names of the types of the metacommunity being analysed

# Arguments:

- `dm`: DiversityMeasure

# Returns:

- Vector of Strings of names of types.
"""
function gettypenames(dm::DiversityMeasure)
    return dm.typenames
end

"""
    getASCIIName(dm::DiversityMeasure)

Return the ASCII name of the DiversityMeasure

# Arguments:

- `dm`: DiversityMeasure

# Returns:

- String containing simple ASCII name of DiversityMeasure
"""
function getASCIIName(dm::DiversityMeasure)
    s = replace(string(typeof(dm)), "Diversity.", "")
    replace(s, r"{.*}$", "")
end

"""
    getName(dm::DiversityMeasure)

Return the character corresponding to the DiversityMeasure.

# Arguments:

- `dm`: DiversityMeasure

# Returns:

- String containing unicode (greek) name of DiversityMeasure.
"""
function getName end

"""
    getFullName(dm::DiversityMeasure)

Return the full name of the DiversityMeasure.

# Arguments:

- `dm`: DiversityMeasure

# Returns:

- String containing full descriptive name of DiversityMeasure
"""
function getFullName end

@compat (dl::DiversityLevel)(dm::DiversityMeasure) = getPartitionFunction(dm, dl)
@compat (dl::DiversityLevel)(dm::DiversityMeasure, qs) = getPartitionFunction(dm, dl)(qs)

"""
    PowerMeanMeasure

This abstract DiversityMeasure subtype is the supertype of all
diversity measures which are straight power means. PowerMeanMeasure
subtypes allow you to calculate and cache any kind of diversity of a
metacommunity.
"""
@compat abstract type PowerMeanMeasure{FP, AbArray} <: DiversityMeasure{FP, AbArray} end

"""
    RelativeEntropyMeasure

This abstract DiversityMeasure subtype is the supertype of all
diversity measures which are relative entropy-based diversity measures.
RelativeEntropyMeasure subtypes allow you to calculate and cache any
kind of diversity of a metacommunity.
"""
@compat abstract type RelativeEntropyMeasure{FP, AbArray} <: DiversityMeasure{FP, AbArray} end

"""
    inddiv(measure::DiversityMeasure, q::Real)
    inddiv(measure::DiversityMeasure, qs::AbstractVector{Real})

Takes a diversity measure and single order or vector of orders, and
returns a DataFrame containing the individual diversities for those values.

# Arguments:

- `dm`: DiversityMeasure
- `q` / `qs`: a single order or a vector of orders

# Returns:

- Returns individual diversities of `dm` for a single order `q` or a
  vector of order `qs`.
"""
function inddiv end

@inline function inddiv(measure::DiversityMeasure, q::Real)
    raw = inddiv_raw(measure, q)
    types = gettypenames(measure)
    scn = getscnames(measure)
    scs = reshape(scn, 1, length(scn))
    dfs = broadcast((div, tn, pn) -> DataFrame(measure=getASCIIName(measure), q=q,
                                               type_level="type", type_name=tn,
                                               partition_level="subcommunity",
                                               partition_name=pn,
                                               diversity=div),
                    raw, types, scs)
    return reduce(append!, dfs)
end

@inline function inddiv(measure::DiversityMeasure, qs::AbstractVector)
    mapreduce(q -> inddiv(measure, q), append!, qs)
end

@inline function inddiv(meta::AbstractMetacommunity, qs)
    mapreduce(dm -> inddiv(dm(meta), qs),
              append!,
              [RawAlpha, NormalisedAlpha,
               RawBeta, NormalisedBeta,
               RawRho, NormalisedRho, Gamma])
end

@inline function inddiv_raw(measure::DiversityMeasure, ::Real)
    measure.diversities
end

"""
    subdiv(measure::DiversityMeasure, q::Real)
    subdiv(measure::DiversityMeasure, qs::AbstractVector{Real})

Takes a diversity measure and single order or vector of orders, and
calculates and returns the subcommunity diversities for those values.

# Arguments:
- `dm`: DiversityMeasure
- `q` / `qs`: a single order or a vector of orders

# Returns:

- Returns subcommunity diversities of `dm` for a single order `q` or a
  vector of order `qs`.
"""
function subdiv end

@inline function subdiv(measure::DiversityMeasure, q::Real)
    raw = subdiv_raw(measure, q)
    scs = getscnames(measure)
    dfs = broadcast((div, pn) -> DataFrame(measure=getASCIIName(measure), q=q,
                                           type_level="types", type_name="",
                                           partition_level="subcommunity",
                                           partition_name=pn,
                                           diversity=div),
                    raw, scs)
    return reduce(append!, dfs)
end

@inline function subdiv(measure::DiversityMeasure, qs::AbstractVector)
    mapreduce(q -> subdiv(measure, q), append!, qs)
end

@inline function subdiv(meta::AbstractMetacommunity, qs)
    mapreduce(dm -> subdiv(dm(meta), qs),
              append!,
              [RawAlpha, NormalisedAlpha,
               RawBeta, NormalisedBeta,
               RawRho, NormalisedRho, Gamma])
end

@inline function subdiv_raw(measure::PowerMeanMeasure, q::Real)
    powermean(inddiv_raw(measure, q), 1.0 - q, measure.abundances)
end

@inline function subdiv_raw(measure::RelativeEntropyMeasure, q::Real)
    powermean(inddiv_raw(measure, q), q - 1.0, measure.abundances)
end

"""
    metadiv(measure::DiversityMeasure, q::Real)
    metadiv(measure::DiversityMeasure, qs::AbstractVector{Real})

Takes a diversity measure and single order or vector of orders, and
calculates and returns the metacommunity diversities for those values.

# Arguments:

- `dm`: DiversityMeasure
- `q` / `qs`: a single order or a vector of orders

# Returns:

- Returns metacommunity diversities of `dm` for a single order `q` or a
  vector of order `qs`.
"""
function metadiv end

@inline function metadiv(measure::DiversityMeasure, q::Real)
    raw = metadiv_raw(measure, q)
    return DataFrame(measure=getASCIIName(measure), q=q,
                     type_level="types", type_name="",
                     partition_level="metacommunity",
                     partition_name="",
                     diversity=raw)
end

@inline function metadiv(measure::DiversityMeasure, qs::AbstractVector)
    mapreduce(q -> metadiv(measure, q), append!, qs)
end

@inline function metadiv(meta::AbstractMetacommunity, qs)
    mapreduce(dm -> metadiv(dm(meta), qs),
              append!,
              [RawAlpha, NormalisedAlpha,
               RawBeta, NormalisedBeta,
               RawRho, NormalisedRho, Gamma])
end

@inline function metadiv_raw(measure::DiversityMeasure, q::Real)
    powermean(vectorise(subdiv_raw(measure, q)), 1.0 - q, measure.weights)
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
    typenames::Vector{String}
    scnames::Vector{String}
end

function RawAlpha(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    w = vectorise(getweight(meta))
    types = gettypes(meta)
    part = getpartition(meta)
    RawAlpha{eltype(ab), typeof(ab)}(ab, w, getordinariness!(meta) .^ -1,
                                     getnames(types), getnames(part))
    
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
    typenames::Vector{String}
    scnames::Vector{String}
end

function NormalisedAlpha(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    ws = getweight(meta)
    w = vectorise(ws)
    types = gettypes(meta)
    part = getpartition(meta)
    NormalisedAlpha{eltype(ab), typeof(ab)}(ab, w,
                                            ws ./ getordinariness!(meta),
                                            getnames(types), getnames(part))
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
    typenames::Vector{String}
    scnames::Vector{String}
end

function RawBeta(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    w = vectorise(getweight(meta))
    types = gettypes(meta)
    part = getpartition(meta)
    RawBeta{eltype(ab), typeof(ab)}(ab, w,
                                    getordinariness!(meta) ./
                                    getmetaordinariness!(meta),
                                    getnames(types), getnames(part))
end

const Distinctiveness = RawBeta

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
    typenames::Vector{String}
    scnames::Vector{String}
end

function NormalisedBeta(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    ws = getweight(meta)
    w = vectorise(ws)
    types = gettypes(meta)
    part = getpartition(meta)
    NormalisedBeta{eltype(ab), typeof(ab)}(ab, w,
                                           getordinariness!(meta) ./
                                           (getmetaordinariness!(meta) .* ws),
                                           getnames(types), getnames(part))
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
    typenames::Vector{String}
    scnames::Vector{String}
end

function RawRho(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    w = vectorise(getweight(meta))
    types = gettypes(meta)
    part = getpartition(meta)
    RawRho{eltype(ab), typeof(ab)}(ab, w,
                                   getmetaordinariness!(meta) ./
                                   getordinariness!(meta),
                                   getnames(types), getnames(part))
end

const Redundancy = RawRho

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
    typenames::Vector{String}
    scnames::Vector{String}
end

function NormalisedRho(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    ws = getweight(meta)
    w = vectorise(ws)
    types = gettypes(meta)
    part = getpartition(meta)
    NormalisedRho{eltype(ab), typeof(ab)}(ab, w,
                                          (getmetaordinariness!(meta) .* ws) ./
                                          getordinariness!(meta),
                                          getnames(types), getnames(part))
end

const Representativeness = NormalisedRho

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
    typenames::Vector{String}
    scnames::Vector{String}
end

function Gamma(meta::AbstractMetacommunity)
    ab = getabundance(meta)
    ws = getweight(meta)
    w = vectorise(ws)
    types = gettypes(meta)
    part = getpartition(meta)
    Gamma{eltype(ab), typeof(ab)}(ab, w,
                                  ones(eltype(ws), size(ws)) ./
                                  getmetaordinariness!(meta),
                                  getnames(types), getnames(part))
end

getName(::Gamma) = "γ"
getFullName(::Gamma) = "gamma diversity"
