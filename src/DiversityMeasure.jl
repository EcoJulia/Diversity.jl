using DataFrames
using EcoBase: AbstractAssemblage

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
abstract type DiversityMeasure{FP <: AbstractFloat,
                               AbMatrix <: AbstractMatrix,
                               DivArray <: AbstractArray,
                               MC <: AbstractAssemblage} end

"""
    getASCIIName(dm::DiversityMeasure)

Return the ASCII name of the DiversityMeasure

# Arguments:

- `dm`: DiversityMeasure

# Returns:

- String containing simple ASCII name of DiversityMeasure
"""
function getASCIIName(dm::DiversityMeasure)
    s = replace(string(typeof(dm)), "Diversity." => "")
    replace(s, r"{.*}$" => "")
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

"""
    _getmeta(dm::DiversityMeasure)

Return the metacommunity belonging to the DiversityMeasure.
"""
function _getmeta end

getsubcommunitynames(dm::DiversityMeasure) = getsubcommunitynames(_getmeta(dm))
gettypenames(dm::DiversityMeasure) = gettypenames(_getmeta(dm))
getdiversityname(dm::DiversityMeasure) = getdiversityname(_getmeta(dm))

(dl::DiversityLevel)(dm::DiversityMeasure) = getPartitionFunction(dm, dl)
(dl::DiversityLevel)(dm::DiversityMeasure, qs) = getPartitionFunction(dm, dl)(qs)

"""
    PowerMeanMeasure

This abstract DiversityMeasure subtype is the supertype of all
diversity measures which are straight power means. PowerMeanMeasure
subtypes allow you to calculate and cache any kind of diversity of a
metacommunity.
"""
abstract type PowerMeanMeasure{FP, AbMatrix, DivArray, MC} <:
    DiversityMeasure{FP, AbMatrix, DivArray, MC} end

"""
    RelativeEntropyMeasure

This abstract DiversityMeasure subtype is the supertype of all
diversity measures which are relative entropy-based diversity measures.
RelativeEntropyMeasure subtypes allow you to calculate and cache any
kind of diversity of a metacommunity.
"""
abstract type RelativeEntropyMeasure{FP, AbMatrix, DivArray, MC} <:
    DiversityMeasure{FP, AbMatrix, DivArray, MC} end

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
    scn = getsubcommunitynames(measure)
    scs = reshape(scn, 1, length(scn))
    dfs = broadcast((div, tn, pn) ->
                    DataFrame(div_type=getdiversityname(measure),
                              measure=getASCIIName(measure),
                              q=q,
                              type_level="type", type_name=tn,
                              partition_level="subcommunity",
                              partition_name=pn,
                              diversity=div),
                    raw, types, scs)
    df = reduce(append!, dfs)
    cols = addedoutputcols(_getmeta(measure))
    if length(cols) > 0
        len = length(df)
        data = getaddedoutput(_getmeta(measure))
        for col in keys(cols)
            insert!(df, length(df) + 1, data[col], col)
        end
    end
    return df
end

@inline function inddiv(measure::DiversityMeasure, qs::AbstractVector)
    mapreduce(q -> inddiv(measure, q), append!, qs)
end

@inline function inddiv(meta::AbstractAssemblage, qs)
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
    scs = getsubcommunitynames(measure)
    dfs = broadcast((div, pn) -> DataFrame(div_type=getdiversityname(measure),
                                           measure=getASCIIName(measure), q=q,
                                           type_level="types", type_name="",
                                           partition_level="subcommunity",
                                           partition_name=pn,
                                           diversity=div),
                    raw, scs)
    df = reduce(append!, dfs)
    cols = addedoutputcols(_getmeta(measure))
    if length(cols) > 0
        len = length(df)
        data = getaddedoutput(_getmeta(measure))
        for col in keys(cols)
            insert!(df, length(df) + 1, data[col], col)
        end
    end
    return df
end

@inline function subdiv(measure::DiversityMeasure, qs::AbstractVector)
    mapreduce(q -> subdiv(measure, q), append!, qs)
end

@inline function subdiv(meta::AbstractAssemblage, qs)
    mapreduce(dm -> subdiv(dm(meta), qs),
              append!,
              [RawAlpha, NormalisedAlpha,
               RawBeta, NormalisedBeta,
               RawRho, NormalisedRho, Gamma])
end

@inline function subdiv_raw(measure::PowerMeanMeasure, q::Real)
    powermean(inddiv_raw(measure, q), one(q) - q, measure.abundances)
end

@inline function subdiv_raw(measure::RelativeEntropyMeasure, q::Real)
    powermean(inddiv_raw(measure, q), q - one(q), measure.abundances)
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
    df = DataFrame(div_type=getdiversityname(measure),
                   measure=getASCIIName(measure), q=q,
                   type_level="types", type_name="",
                   partition_level="metacommunity",
                   partition_name="",
                   diversity=raw)
   cols = addedoutputcols(_getmeta(measure))
   if length(cols) > 0
       len = length(df)
       data = getaddedoutput(_getmeta(measure))
       for col in keys(cols)
           insert!(df, length(df) + 1, data[col], col)
       end
   end
   return df
end

@inline function metadiv(measure::DiversityMeasure, qs::AbstractVector)
    mapreduce(q -> metadiv(measure, q), append!, qs)
end

@inline function metadiv(meta::AbstractAssemblage, qs)
    mapreduce(dm -> metadiv(dm(meta), qs),
              append!,
              [RawAlpha, NormalisedAlpha,
               RawBeta, NormalisedBeta,
               RawRho, NormalisedRho, Gamma])
end

@inline function metadiv_raw(measure::DiversityMeasure, q::Real)
    powermean(subdiv_raw(measure, q), one(q) - q, measure.weights)
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
    RawAlpha

Calculates raw alpha diversity (α) of all of the individuals in a
metacommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
struct RawAlpha{FP, AbMatrix, DivArray, MC} <:
    PowerMeanMeasure{FP, AbMatrix, DivArray, MC}
    abundances::AbMatrix
    weights::Vector{FP}
    diversities::DivArray
    meta::MC
end

function RawAlpha(meta::M) where M <: AbstractAssemblage
    ab = getabundance(meta)
    ws = getweight(meta)
    value = getordinariness!(meta) .^ -1
    return RawAlpha{eltype(ab), typeof(ab),
                    typeof(value), M}(ab, ws, value, meta)
end

getName(::RawAlpha) = "α"
getFullName(::RawAlpha) = "raw alpha diversity"
_getmeta(m::RawAlpha) = m.meta

"""
    NormalisedAlpha

Calculates normalised alpha diversity (ᾱ) of all of the individuals in
a metacommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
struct NormalisedAlpha{FP, AbMatrix, DivArray, MC} <:
    PowerMeanMeasure{FP, AbMatrix, DivArray, MC}
    abundances::AbMatrix
    weights::Vector{FP}
    diversities::DivArray
    meta::MC
end

function NormalisedAlpha(meta::M) where M <: AbstractAssemblage
    ab = getabundance(meta)
    ws = getweight(meta)
    value = ws' ./ getordinariness!(meta)
    return NormalisedAlpha{eltype(ab), typeof(ab),
                           typeof(value), M}(ab, ws, value, meta)
end

getName(::NormalisedAlpha) = "ᾱ"
getFullName(::NormalisedAlpha) = "normalised alpha diversity"
_getmeta(m::NormalisedAlpha) = m.meta

"""
    RawBeta

Calculates distinctiveness (β, raw beta diversity) of all of the individuals in a
metacommunity, and caches them for subsequent analysis. This is a
subtype of RelativeEntropyMeasure, meaning that subcommunity and type
composite diversity measures are relative entropies, and their
composite types are powermeans of those measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
struct RawBeta{FP, AbMatrix, DivArray, MC} <:
    RelativeEntropyMeasure{FP, AbMatrix, DivArray, MC}
    abundances::AbMatrix
    weights::Vector{FP}
    diversities::DivArray
    meta::MC
end

function RawBeta(meta::M) where M <: AbstractAssemblage
    ab = getabundance(meta)
    ws = getweight(meta)
    value = getordinariness!(meta) ./ getmetaordinariness!(meta)
    return RawBeta{eltype(ab), typeof(ab),
                   typeof(value), M}(ab, ws, value, meta)
end

const Distinctiveness = RawBeta

getName(::RawBeta) = "β"
getFullName(::RawBeta) = "distinctiveness"
_getmeta(m::RawBeta) = m.meta

"""
    NormalisedBeta

Calculates normalised beta diversity (β̄) of all of the individuals in
a metacommunity, and caches them for subsequent analysis. This is a
subtype of RelativeEntropyMeasure, meaning that subcommunity and type
composite diversity measures are relative entropies, and their
composite types are powermeans of those measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
struct NormalisedBeta{FP, AbMatrix, DivArray, MC} <:
    RelativeEntropyMeasure{FP, AbMatrix, DivArray, MC}
    abundances::AbMatrix
    weights::Vector{FP}
    diversities::DivArray
    meta::MC
end

function NormalisedBeta(meta::M) where M <: AbstractAssemblage
    ab = getabundance(meta)
    ws = getweight(meta)
    value = getordinariness!(meta) ./ (getmetaordinariness!(meta) .* ws')
    return NormalisedBeta{eltype(ab), typeof(ab),
                          typeof(value), M}(ab, ws, value, meta)
end

getName(::NormalisedBeta) = "β̄"
getFullName(::NormalisedBeta) = "effective number of subcommunities"
_getmeta(m::NormalisedBeta) = m.meta

"""
    RawRho

Calculates redundancy (ρ, raw beta diversity) of all of the
individuals in a metacommunity, and caches them for subsequent
analysis. This is a subtype of PowerMeanMeasure, meaning that all
composite diversity measures are simple powermeans of the individual
measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
struct RawRho{FP, AbMatrix, DivArray, MC} <:
    PowerMeanMeasure{FP, AbMatrix, DivArray, MC}
    abundances::AbMatrix
    weights::Vector{FP}
    diversities::DivArray
    meta::MC
end

function RawRho(meta::M) where M <: AbstractAssemblage
    ab = getabundance(meta)
    ws = getweight(meta)
    value = getmetaordinariness!(meta) ./ getordinariness!(meta)
    return RawRho{eltype(ab), typeof(ab),
                  typeof(value), M}(ab, ws, value, meta)
end

const Redundancy = RawRho

getName(::RawRho) = "ρ"
getFullName(::RawRho) = "redundancy"
_getmeta(m::RawRho) = m.meta

"""
    NormalisedRho

Calculates redundancy (ρ̄, normalised beta diversity) of all of the
individuals in a metacommunity, and caches them for subsequent
analysis. This is a subtype of PowerMeanMeasure, meaning that all
composite diversity measures are simple powermeans of the individual
measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
struct NormalisedRho{FP, AbMatrix, DivArray, MC} <:
    PowerMeanMeasure{FP, AbMatrix, DivArray, MC}
    abundances::AbMatrix
    weights::Vector{FP}
    diversities::DivArray
    meta::MC
end

function NormalisedRho(meta::M) where M <: AbstractAssemblage
    ab = getabundance(meta)
    ws = getweight(meta)
    value = (getmetaordinariness!(meta) .* ws') ./ getordinariness!(meta)
    return NormalisedRho{eltype(ab), typeof(ab),
                         typeof(value), M}(ab, ws, value, meta)
end

const Representativeness = NormalisedRho

getName(::NormalisedRho) = "ρ̄"
getFullName(::NormalisedRho) = "representativeness"
_getmeta(m::NormalisedRho) = m.meta

"""
    Gamma

Calculates gamma diversity (γ) of all of the individuals in a
metacommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

#### Constructor arguments:

- `meta`: a Metacommunity
"""
struct Gamma{FP, AbMatrix, DivArray, MC} <:
    PowerMeanMeasure{FP, AbMatrix, DivArray, MC}
    abundances::AbMatrix
    weights::Vector{FP}
    diversities::DivArray
    meta::MC
end

function Gamma(meta::M) where M <: AbstractAssemblage
    ab = getabundance(meta)
    ws = getweight(meta)
    value = fill!(similar(ws), 1)' ./ getmetaordinariness!(meta)
    return Gamma{eltype(ab), typeof(ab), typeof(value), M}(ab, ws, value, meta)
end

getName(::Gamma) = "γ"
getFullName(::Gamma) = "gamma diversity"
_getmeta(m::Gamma) = m.meta

RecipesBase.@recipe function f(var::Tuple{<: DiversityMeasure,
                                          <: Real})
    title := getFullName(var[1]) * " (q = $(var[2])) - " * getdiversityname(var[1]) *
        " diversity"
    colorbar_title := getASCIIName(var[1])
    subdiv(var...)[:diversity], getcoords(places(_getmeta(var[1])))
end
