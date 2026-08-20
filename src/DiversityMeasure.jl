# SPDX-License-Identifier: BSD-2-Clause

using DataFrames
using Tables
import EcoBase: AbstractAssemblage

"""
### Enumeration of levels that can exist / be calculated for a metacommunity.
"""
@enum DiversityLevel individualDiversity subcommunityDiversity communityDiversity typeDiversity typeCollectionDiversity metacommunityDiversity

# The seven measures, in the order every "all measures" call has always produced them. A function
# rather than a const because the structs are defined further down this file.
function _allmeasures()
    return (RawAlpha, NormalisedAlpha, RawBeta, NormalisedBeta,
            RawRho, NormalisedRho, Gamma)
end

# Seven of a result's eight columns are repeated patterns - only `diversity` is unpredictable data.
#
# Held as rules instead they cost nothing, and nothing downstream has to know. `Tables.materializer`
# decides: `DataFrame` copies its columns by default, and `copy` of an `AbstractVector` goes through
# `similar` and `copyto!`, so a caller asking for a DataFrame gets ordinary `Vector`s exactly as
# before.
#
# Warning: `_vcatcolumns` concatenates with `vcat`, which materialises, so a result covering several
# orders, measures or levels loses this again. Streaming those as partitions is a separate question.

# One value, repeated for every row.
struct ConstantColumn{T} <: AbstractVector{T}
    value::T
    len::Int
end

Base.size(col::ConstantColumn) = (col.len,)
Base.IndexStyle(::Type{<:ConstantColumn}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(col::ConstantColumn, ::Int) = col.value

# A short list cycled to fill a column, replacing `repeat`. `runlength` is how many consecutive rows
# share an entry before the next is used -- 1 for the type names, which cycle fastest, and ntypes
# for the subcommunity names, which change once per block -- so this covers `repeat(v, outer = k)`
# and `repeat(v, inner = k)` alike.
#
# The inner constructor copies, so the column can never alias the names it was built from. That is
# free relative to what it replaces: the list is ntypes or nplaces long where the column is their
# product.
struct RepeatedColumn{T, V <: AbstractVector{T}} <: AbstractVector{T}
    values::V
    runlength::Int
    len::Int

    function RepeatedColumn(values::V, runlength::Int,
                            len::Int) where {T, V <: AbstractVector{T}}
        return new{T, V}(copy(values), runlength, len)
    end
end

Base.size(col::RepeatedColumn) = (col.len,)
Base.IndexStyle(::Type{<:RepeatedColumn}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(col::RepeatedColumn, i::Int)
    return col.values[mod1(cld(i, col.runlength), length(col.values))]
end

# The columns of a result, as a NamedTuple of equal-length vectors. That is already a Tables source,
# so it can be handed to `Tables.materializer(sink)` for any table type the caller asks for -- a
# DataFrame by default, but equally an Arrow table, a CSV sink or anything else implementing the
# interface. Building the columns whole rather than a row at a time is also what makes this cheap;
# see the performance notes in CLAUDE.md.
#
# `addedoutputcols` lets a types object contribute extra columns (the Phylo extension adds
# `:treename`), merged in here rather than inserted afterwards, since a NamedTuple is immutable and
# a sink may not support insertion at all.
function _addedcolumns(measure, columns, n)
    cols = addedoutputcols(_getmeta(measure))
    isempty(cols) && return columns
    data = getaddedoutput(_getmeta(measure))
    extra = NamedTuple(col => ConstantColumn(data[col], n)
                       for col in keys(cols))
    return merge(columns, extra)
end

# Which set of columns a DiversityLevel asks for. The counterpart of `getPartitionFunction`, but
# returning the columns rather than a materialised table, so that a caller wanting several levels at
# once builds only one.
function _levelcolumns(level::DiversityLevel, measure, qs)
    level == individualDiversity && return _inddiv_columns(measure, qs)
    level == subcommunityDiversity && return _subdiv_columns(measure, qs)
    level == metacommunityDiversity && return _metadiv_columns(measure, qs)
    return error("Can't calculate diversity for $level")
end

# Concatenate several results' columns, which is what asking for several orders, measures or levels
# at once produces. Done on the columns so that only one table is ever materialised.
function _vcatcolumns(parts)
    length(parts) == 1 && return only(parts)
    ks = keys(first(parts))
    return NamedTuple{ks}(map(k -> reduce(vcat, (part[k] for part in parts)),
                              ks))
end

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

# The individual diversities of a measure -- one number for every type in every subcommunity --
# held as a rule for computing an element rather than as an array of them. Every one of the seven
# measures is an elementwise combination of at most three things: the ordinariness of the
# individuals of a type in a subcommunity, the ordinariness of that type in the metacommunity, and
# the subcommunity's weight. Only the first is large, and the measure already has it -- the
# metacommunity caches it -- so materialising the combination doubles the memory for no new
# information.
#
# This is a trade, not a free win, and the direction is worth being plain about. Computing an
# element reads exactly the same memory the array would have -- the ordinariness, either way -- but
# adds one division, and saves storing the whole array.
#
# The rule is a closure so that each measure captures exactly the arrays it uses -- alpha and gamma
# do not need both ordinarinesses, and forcing them to would compute one they have no use for.
struct IndividualDiversities{FP <: AbstractFloat, F} <: AbstractMatrix{FP}
    value::F
    dims::Tuple{Int, Int}

    function IndividualDiversities{FP}(value::F,
                                       dims::Tuple{Int, Int}) where
        {FP <: AbstractFloat, F}
        return new{FP, F}(value, dims)
    end
end

Base.size(divs::IndividualDiversities) = divs.dims
Base.IndexStyle(::Type{<:IndividualDiversities}) = IndexCartesian()
Base.@propagate_inbounds function Base.getindex(divs::IndividualDiversities,
                                                i::Int, j::Int)
    return divs.value(i, j)
end

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
    return string(nameof(typeof(dm)))
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
function _getmeta(dm::DiversityMeasure)
    return dm.meta
end

getsubcommunitynames(dm::DiversityMeasure) = getsubcommunitynames(_getmeta(dm))
gettypenames(dm::DiversityMeasure) = gettypenames(_getmeta(dm))
getdiversityname(dm::DiversityMeasure) = getdiversityname(_getmeta(dm))

(dl::DiversityLevel)(dm::DiversityMeasure) = getPartitionFunction(dm, dl)
function (dl::DiversityLevel)(dm::DiversityMeasure, qs)
    return getPartitionFunction(dm, dl)(qs)
end

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

function _inddiv_columns(measure::DiversityMeasure, q::Real)
    raw = inddiv_raw(measure, q)
    types = gettypenames(measure)
    scn = getsubcommunitynames(measure)
    nt, ns = length(types), length(scn)
    n = nt * ns
    # The individual diversities are held as a rule rather than an array, so this is where they
    # are materialised -- which has to happen anyway, since a column of the output is a vector.
    divs = Matrix{eltype(raw)}(undef, nt, ns)
    divs .= raw
    # The row order `reduce(append!, ...)` over a column-major matrix used to produce: types
    # cycling fastest within each subcommunity.
    columns = (div_type = ConstantColumn(getdiversityname(measure), n),
               measure = ConstantColumn(getASCIIName(measure), n),
               q = ConstantColumn(q, n),
               type_level = ConstantColumn("type", n),
               type_name = RepeatedColumn(types, 1, n),
               partition_level = ConstantColumn("subcommunity", n),
               partition_name = RepeatedColumn(scn, nt, n),
               diversity = vec(divs))
    return _addedcolumns(measure, columns, n)
end

function _inddiv_columns(measure::DiversityMeasure, qs::AbstractVector)
    return _vcatcolumns([_inddiv_columns(measure, q) for q in qs])
end

function _inddiv_columns(meta::AbstractAssemblage, qs)
    return _vcatcolumns([_inddiv_columns(dm(meta), qs) for dm in _allmeasures()])
end

@inline function inddiv(sink, measure::DiversityMeasure, qs)
    return Tables.materializer(sink)(_inddiv_columns(measure, qs))
end

@inline function inddiv(sink, meta::AbstractAssemblage, qs)
    return Tables.materializer(sink)(_inddiv_columns(meta, qs))
end

@inline inddiv(measure::DiversityMeasure, qs) = inddiv(DataFrame, measure, qs)
@inline inddiv(meta::AbstractAssemblage, qs) = inddiv(DataFrame, meta, qs)

@inline function inddiv_raw(measure::DiversityMeasure, ::Real)
    return measure.diversities
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

function _subdiv_columns(measure::DiversityMeasure, q::Real)
    raw = subdiv_raw(measure, q)
    scn = getsubcommunitynames(measure)
    n = length(scn)
    divs = Vector{eltype(raw)}(undef, n)
    divs .= raw
    columns = (div_type = ConstantColumn(getdiversityname(measure), n),
               measure = ConstantColumn(getASCIIName(measure), n),
               q = ConstantColumn(q, n),
               type_level = ConstantColumn("types", n),
               type_name = ConstantColumn("", n),
               partition_level = ConstantColumn("subcommunity", n),
               partition_name = copy(scn),
               diversity = divs)
    return _addedcolumns(measure, columns, n)
end

function _subdiv_columns(measure::DiversityMeasure, qs::AbstractVector)
    return _vcatcolumns([_subdiv_columns(measure, q) for q in qs])
end

function _subdiv_columns(meta::AbstractAssemblage, qs)
    return _vcatcolumns([_subdiv_columns(dm(meta), qs) for dm in _allmeasures()])
end

@inline function subdiv(sink, measure::DiversityMeasure, qs)
    return Tables.materializer(sink)(_subdiv_columns(measure, qs))
end

@inline function subdiv(sink, meta::AbstractAssemblage, qs)
    return Tables.materializer(sink)(_subdiv_columns(meta, qs))
end

@inline subdiv(measure::DiversityMeasure, qs) = subdiv(DataFrame, measure, qs)
@inline subdiv(meta::AbstractAssemblage, qs) = subdiv(DataFrame, meta, qs)

@inline function subdiv_raw(measure::PowerMeanMeasure, q::Real)
    return powermean(inddiv_raw(measure, q), one(q) - q, measure.abundances)
end

@inline function subdiv_raw(measure::RelativeEntropyMeasure, q::Real)
    return powermean(inddiv_raw(measure, q), q - one(q), measure.abundances)
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

function _metadiv_columns(measure::DiversityMeasure, q::Real)
    raw = metadiv_raw(measure, q)
    columns = (div_type = [getdiversityname(measure)],
               measure = [getASCIIName(measure)],
               q = [q],
               type_level = ["types"],
               type_name = [""],
               partition_level = ["metacommunity"],
               partition_name = [""],
               diversity = [raw])
    return _addedcolumns(measure, columns, 1)
end

function _metadiv_columns(measure::DiversityMeasure, qs::AbstractVector)
    return _vcatcolumns([_metadiv_columns(measure, q) for q in qs])
end

function _metadiv_columns(meta::AbstractAssemblage, qs)
    return _vcatcolumns([_metadiv_columns(dm(meta), qs)
                         for dm in _allmeasures()])
end

@inline function metadiv(sink, measure::DiversityMeasure, qs)
    return Tables.materializer(sink)(_metadiv_columns(measure, qs))
end

@inline function metadiv(sink, meta::AbstractAssemblage, qs)
    return Tables.materializer(sink)(_metadiv_columns(meta, qs))
end

@inline metadiv(measure::DiversityMeasure, qs) = metadiv(DataFrame, measure, qs)
@inline metadiv(meta::AbstractAssemblage, qs) = metadiv(DataFrame, meta, qs)

@inline function metadiv_raw(measure::DiversityMeasure, q::Real)
    return powermean(subdiv_raw(measure, q), one(q) - q, measure.weights)
end

function getPartitionFunction(measure::DiversityMeasure, level::DiversityLevel)
    if (level == individualDiversity)
        return function (qs)
            return inddiv(measure, qs)
        end
    elseif (level == subcommunityDiversity)
        function (qs)
            return subdiv(measure, qs)
        end
    elseif (level == metacommunityDiversity)
        function (qs)
            return metadiv(measure, qs)
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

Per subcommunity, it is an estimate of naive-community metacommunity diversity
— the diversity the whole metacommunity would have if this subcommunity shared
no types, and no similarity, with any other. Averaged over the subcommunities it
gives naive-community metacommunity diversity itself, which is an upper bound on
the true metacommunity diversity `Gamma`. It is `NormalisedAlpha` measured per
individual rather than per subcommunity.

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

function RawAlpha(meta::M) where {M <: AbstractAssemblage}
    ab = getabundance(meta)
    ws = getweight(meta)
    zp = getordinariness!(meta)
    value = IndividualDiversities{eltype(ab)}((i, j) -> inv(zp[i, j]),
                                              size(ab))
    return RawAlpha{eltype(ab), typeof(ab),
                    typeof(value), M}(ab, ws, value, meta)
end

getName(::RawAlpha) = "α"
getFullName(::RawAlpha) = "estimate of naive-community metacommunity diversity"

"""
    NormalisedAlpha

Calculates normalised alpha diversity (ᾱ) of all of the individuals in
a metacommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

Per subcommunity, it is the similarity-sensitive diversity of that subcommunity
in isolation — what its diversity would be if it were the whole community.
Averaged over the subcommunities it gives their average diversity, which is
invariant under shattering.

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

function NormalisedAlpha(meta::M) where {M <: AbstractAssemblage}
    ab = getabundance(meta)
    ws = getweight(meta)
    zp = getordinariness!(meta)
    value = IndividualDiversities{eltype(ab)}((i, j) -> ws[j] / zp[i, j],
                                              size(ab))
    return NormalisedAlpha{eltype(ab), typeof(ab),
                           typeof(value), M}(ab, ws, value, meta)
end

getName(::NormalisedAlpha) = "ᾱ"
getFullName(::NormalisedAlpha) = "diversity of subcommunity in isolation"

"""
    RawBeta

Calculates distinctiveness (β, raw beta diversity) of all of the individuals in a
metacommunity, and caches them for subsequent analysis. This is a
subtype of RelativeEntropyMeasure, meaning that subcommunity and type
composite diversity measures are relative entropies, and their
composite types are powermeans of those measures.

Per subcommunity, it is the **distinctiveness** of that subcommunity: how much
of it is unlike the rest of the metacommunity, whether through types found
nowhere else or through low similarity to the types that are. It reaches its
maximum of 1 when every individual in the subcommunity is completely dissimilar
to every individual outside it, and is small when the subcommunity has much in
common with the rest. Averaged over the subcommunities it gives their average
distinctiveness, which can be read as a kind of turnover. It is the reciprocal
of `RawRho`.

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

function RawBeta(meta::M) where {M <: AbstractAssemblage}
    ab = getabundance(meta)
    ws = getweight(meta)
    zp = getordinariness!(meta)
    zP = getmetaordinariness!(meta)
    value = IndividualDiversities{eltype(ab)}((i, j) -> zp[i, j] / zP[i],
                                              size(ab))
    return RawBeta{eltype(ab), typeof(ab),
                   typeof(value), M}(ab, ws, value, meta)
end

# The paper's own name for this measure; provided so that code can read the way the framework
# describes it. Not exported - reach it as `Diversity.Distinctiveness`.
const Distinctiveness = RawBeta

getName(::RawBeta) = "β"
getFullName(::RawBeta) = "distinctiveness"

"""
    NormalisedBeta

Calculates normalised beta diversity (β̄) of all of the individuals in
a metacommunity, and caches them for subsequent analysis. This is a
subtype of RelativeEntropyMeasure, meaning that subcommunity and type
composite diversity measures are relative entropies, and their
composite types are powermeans of those measures.

Per subcommunity, it is an estimate of the effective number of distinct
subcommunities, and is high when a subcommunity is both distinctive and small.
Averaged over the subcommunities it gives the effective number of distinct
subcommunities itself, which is at most the number of subcommunities — reaching
that maximum when they are completely distinct and of equal size — and which is
invariant under shattering. It is the reciprocal of `NormalisedRho`.

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

function NormalisedBeta(meta::M) where {M <: AbstractAssemblage}
    ab = getabundance(meta)
    ws = getweight(meta)
    zp = getordinariness!(meta)
    zP = getmetaordinariness!(meta)
    value = IndividualDiversities{eltype(ab)}((i, j) -> zp[i, j] /
                                                        (zP[i] * ws[j]),
                                              size(ab))
    return NormalisedBeta{eltype(ab), typeof(ab),
                          typeof(value), M}(ab, ws, value, meta)
end

getName(::NormalisedBeta) = "β̄"
function getFullName(::NormalisedBeta)
    return "estimate of effective number of distinct subcommunities"
end

"""
    RawRho

Calculates redundancy (ρ, raw beta diversity) of all of the
individuals in a metacommunity, and caches them for subsequent
analysis. This is a subtype of PowerMeanMeasure, meaning that all
composite diversity measures are simple powermeans of the individual
measures.

Per subcommunity, it is the **redundancy** of that subcommunity: the extent to
which the diversity of the metacommunity would be preserved if the subcommunity
were lost. It takes its minimum of 1 when nothing resembling the subcommunity
remains elsewhere, so that losing it would lose its diversity entirely. Averaged
over the subcommunities it gives their average redundancy, which rises towards
the *effective* number of subcommunities — the Hill number of their weights — as
they become more alike, reaching the number of subcommunities itself only when
they are also of equal size. It is the reciprocal of `RawBeta`.

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

function RawRho(meta::M) where {M <: AbstractAssemblage}
    ab = getabundance(meta)
    ws = getweight(meta)
    zp = getordinariness!(meta)
    zP = getmetaordinariness!(meta)
    value = IndividualDiversities{eltype(ab)}((i, j) -> zP[i] / zp[i, j],
                                              size(ab))
    return RawRho{eltype(ab), typeof(ab),
                  typeof(value), M}(ab, ws, value, meta)
end

# The paper's own name for this measure. Not exported - reach it as `Diversity.Redundancy`.
const Redundancy = RawRho

getName(::RawRho) = "ρ"
getFullName(::RawRho) = "redundancy"

"""
    NormalisedRho

Calculates representativeness (ρ̄, normalised beta diversity) of all of the
individuals in a metacommunity, and caches them for subsequent
analysis. This is a subtype of PowerMeanMeasure, meaning that all
composite diversity measures are simple powermeans of the individual
measures.

Per subcommunity, it is the **representativeness** of that subcommunity: how
typical it is of the metacommunity as a whole. Where all types are equally
abundant, a subcommunity holding a fraction `r` of them has representativeness
exactly `r` — whatever fraction of the *individuals* it holds, since being the
normalised measure it has the subcommunity's weight divided out. Averaged over
the subcommunities it gives their average
representativeness. In the naive-type case representativeness is at most 1,
attained when the subcommunity has the same type distribution as the
metacommunity — but that bound does **not** hold for a general similarity
matrix. It is the reciprocal of `NormalisedBeta`.

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

function NormalisedRho(meta::M) where {M <: AbstractAssemblage}
    ab = getabundance(meta)
    ws = getweight(meta)
    zp = getordinariness!(meta)
    zP = getmetaordinariness!(meta)
    value = IndividualDiversities{eltype(ab)}((i, j) -> zP[i] * ws[j] /
                                                        zp[i, j],
                                              size(ab))
    return NormalisedRho{eltype(ab), typeof(ab),
                         typeof(value), M}(ab, ws, value, meta)
end

# The paper's own name for this measure. Not exported - reach it as
# `Diversity.Representativeness`.
const Representativeness = NormalisedRho

getName(::NormalisedRho) = "ρ̄"
getFullName(::NormalisedRho) = "representativeness"

"""
    Gamma

Calculates gamma diversity (γ) of all of the individuals in a
metacommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

The two scales read differently here, and the difference matters. Per
subcommunity, it is the contribution *per individual* toward metacommunity
diversity, combining a subcommunity's own diversity with the rarity of its types
in the metacommunity — so a subcommunity of a few very rare types contributes
heavily however dull it looks in isolation. Averaged over the subcommunities it
gives the metacommunity's own similarity-sensitive diversity, the diversity of
the whole taken without regard to how it is divided.

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

function Gamma(meta::M) where {M <: AbstractAssemblage}
    ab = getabundance(meta)
    ws = getweight(meta)
    zP = getmetaordinariness!(meta)
    value = IndividualDiversities{eltype(ab)}((i, j) -> inv(zP[i]), size(ab))
    return Gamma{eltype(ab), typeof(ab), typeof(value), M}(ab, ws, value, meta)
end

getName(::Gamma) = "γ"
function getFullName(::Gamma)
    return "contribution per individual toward metacommunity diversity"
end

RecipesBase.@recipe function f(var::Tuple{<:DiversityMeasure,
                                          <:Real})
    title := getFullName(var[1]) * " (q = $(var[2])) - " *
             getdiversityname(var[1]) *
             " diversity"
    colorbar_title := getASCIIName(var[1])
    return subdiv(var...)[!, :diversity],
           getcoords(places(_getmeta(var[1])))
end
