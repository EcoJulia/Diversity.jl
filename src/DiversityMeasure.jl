# SPDX-License-Identifier: BSD-2-Clause

using DataFrames
import EcoBase: AbstractAssemblage

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

@inline function inddiv(measure::DiversityMeasure, q::Real)
    raw = inddiv_raw(measure, q)
    types = gettypenames(measure)
    scn = getsubcommunitynames(measure)
    nt, ns = length(types), length(scn)
    n = nt * ns
    # Broadcast into the full shape rather than reshaping `raw` directly: most measures hold an
    # ntypes x nsubcommunities array of individual diversities, but Gamma holds one column and
    # relies on it broadcasting across the subcommunities.
    divs = Matrix{eltype(raw)}(undef, nt, ns)
    divs .= raw
    # Built as whole columns, in the order `reduce(append!, ...)` over a column-major matrix used to
    # produce: types cycling fastest within each subcommunity. Every column but the last is either
    # constant or a repetition, so none of them needs to be assembled a row at a time.
    df = DataFrame(div_type = fill(getdiversityname(measure), n),
                   measure = fill(getASCIIName(measure), n),
                   q = fill(q, n),
                   type_level = fill("type", n),
                   type_name = repeat(types, outer = ns),
                   partition_level = fill("subcommunity", n),
                   partition_name = repeat(scn, inner = nt),
                   diversity = vec(divs))
    cols = addedoutputcols(_getmeta(measure))
    if length(cols) > 0
        data = getaddedoutput(_getmeta(measure))
        for col in keys(cols)
            insertcols!(df, ncol(df) + 1, col => data[col])
        end
    end
    return df
end

@inline function inddiv(measure::DiversityMeasure, qs::AbstractVector)
    return mapreduce(q -> inddiv(measure, q), append!, qs)
end

@inline function inddiv(meta::AbstractAssemblage, qs)
    return mapreduce(dm -> inddiv(dm(meta), qs),
                     append!,
                     [RawAlpha, NormalisedAlpha,
                         RawBeta, NormalisedBeta,
                         RawRho, NormalisedRho, Gamma])
end

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

@inline function subdiv(measure::DiversityMeasure, q::Real)
    raw = subdiv_raw(measure, q)
    scn = getsubcommunitynames(measure)
    n = length(scn)
    divs = Vector{eltype(raw)}(undef, n)
    divs .= raw
    df = DataFrame(div_type = fill(getdiversityname(measure), n),
                   measure = fill(getASCIIName(measure), n),
                   q = fill(q, n),
                   type_level = fill("types", n),
                   type_name = fill("", n),
                   partition_level = fill("subcommunity", n),
                   partition_name = copy(scn),
                   diversity = divs)
    cols = addedoutputcols(_getmeta(measure))
    if length(cols) > 0
        data = getaddedoutput(_getmeta(measure))
        for col in keys(cols)
            insertcols!(df, ncol(df) + 1, col => data[col])
        end
    end
    return df
end

@inline function subdiv(measure::DiversityMeasure, qs::AbstractVector)
    return mapreduce(q -> subdiv(measure, q), append!, qs)
end

@inline function subdiv(meta::AbstractAssemblage, qs)
    return mapreduce(dm -> subdiv(dm(meta), qs),
                     append!,
                     [RawAlpha, NormalisedAlpha,
                         RawBeta, NormalisedBeta,
                         RawRho, NormalisedRho, Gamma])
end

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

@inline function metadiv(measure::DiversityMeasure, q::Real)
    raw = metadiv_raw(measure, q)
    df = DataFrame(div_type = getdiversityname(measure),
                   measure = getASCIIName(measure), q = q,
                   type_level = "types", type_name = "",
                   partition_level = "metacommunity",
                   partition_name = "",
                   diversity = raw)
    cols = addedoutputcols(_getmeta(measure))
    if length(cols) > 0
        data = getaddedoutput(_getmeta(measure))
        for col in keys(cols)
            insertcols!(df, ncol(df) + 1, col => data[col])
        end
    end
    return df
end

@inline function metadiv(measure::DiversityMeasure, qs::AbstractVector)
    return mapreduce(q -> metadiv(measure, q), append!, qs)
end

@inline function metadiv(meta::AbstractAssemblage, qs)
    return mapreduce(dm -> metadiv(dm(meta), qs),
                     append!,
                     [RawAlpha, NormalisedAlpha,
                         RawBeta, NormalisedBeta,
                         RawRho, NormalisedRho, Gamma])
end

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
    value = getordinariness!(meta) .^ -1
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
    value = ws' ./ getordinariness!(meta)
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
    value = getordinariness!(meta) ./ getmetaordinariness!(meta)
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
    value = getordinariness!(meta) ./ (getmetaordinariness!(meta) .* ws')
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
    value = getmetaordinariness!(meta) ./ getordinariness!(meta)
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
    value = (getmetaordinariness!(meta) .* ws') ./ getordinariness!(meta)
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
    value = fill!(similar(ws), 1)' ./ getmetaordinariness!(meta)
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
