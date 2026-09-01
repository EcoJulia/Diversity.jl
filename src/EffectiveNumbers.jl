# SPDX-License-Identifier: BSD-2-Clause

"""
    powermean

# Calculates the weighted powermean of a series of numbers

Calculates *order*th power mean of *values*, weighted by
*weights*. By default, *weights* are equal and *order*
is 1, so this is just the arithmetic mean.

# Arguments:

- `values`: values for which to calculate mean
- `order[s]`: order[s] of power mean
- `weights`: weights of elements, normalised to 1 inside function

# Returns:

- weighted power mean(s)
"""
function powermean end
function powermean(values::V1, order::R = 1,
                   weights::V2 = fill!(similar(values), 1)) where
    {R <: Real, FP <: AbstractFloat,
     V1 <: AbstractVector{FP}, V2 <: AbstractVector{FP}}
    length(values) == length(weights) ||
        throw(DimensionMismatch("powermean: Weight and value vectors must be the same length"))

    # Check whether all weights are zero in group.
    # In that case we want to return a NaN
    if iszero(weights)
        return convert(FP, NaN)
    end

    if isinf(order)
        if order > 0 # +Inf -> Maximum
            s = zero(FP)
            for i in eachindex(values, weights)
                @inbounds if (weights[i] > eps(FP)) & (values[i] > s)
                    s = values[i]
                end
            end
            return s
        else # -Inf -> Minimum
            s = convert(FP, Inf)
            for i in eachindex(values, weights)
                @inbounds if (weights[i] > eps(FP)) & (values[i] < s)
                    s = values[i]
                end
            end
            return s
        end
    else
        if order ≈ zero(order)
            s = zero(FP)
            for i in eachindex(values, weights)
                @inbounds if weights[i] > eps(FP)
                    s += weights[i] * log(values[i])
                end
            end
            return exp(s / sum(weights))
        else
            s = zero(FP)
            for i in eachindex(values, weights)
                @inbounds if weights[i] > eps(FP)
                    s += weights[i] * values[i]^order
                end
            end
            return (s / sum(weights))^(one(FP) / order)
        end
    end
end

# This is the next most common case - a vector of orders
function powermean(values::V1,
                   orders::VR,
                   weights::V2 = fill!(similar(values), 1.0)) where
    {R <: Real, VR <: AbstractVector{R},
     FP <: AbstractFloat, V1 <: AbstractVector{FP},
     V2 <: AbstractVector{FP}}
    return map(order -> powermean(values, order, weights), orders)
end

# Whether a subcommunity is empty, given the per-subcommunity weights if the caller had them to
# hand. Without them every column has to find out for itself, by scanning.
_isemptycolumn(::Nothing, ::Int) = false
_isemptycolumn(colweights, col::Int) = iszero(colweights[col])

# The answer for an empty subcommunity, shaped to match what a real column would have returned:
# one number per order.
_emptymean(::Type{FP}, ::Real) where {FP} = convert(FP, NaN)
function _emptymean(::Type{FP}, orders::AbstractVector) where {FP}
    return fill(convert(FP, NaN), length(orders))
end

# This is the next most simple case - matrices with subcommunities, and an order or orders
#
# `colweights` is optional, and is purely an optimisation: it is the weight of each subcommunity,
# which a `Metacommunity` already caches and which the measures therefore have to hand. A
# subcommunity of zero weight holds no individuals, so its diversity is NaN.
function powermean(values::M1, orders,
                   weights::M2 = fill!(similar(values), 1),
                   colweights = nothing) where
    {FP <: AbstractFloat, M1 <: AbstractMatrix{FP},
     M2 <: AbstractMatrix{FP}}
    size(values) == size(weights) ||
        throw(DimensionMismatch("powermean: Weight and value matrixes " *
                                "must be the same size"))
    isnothing(colweights) || length(colweights) == size(values, 2) ||
        throw(DimensionMismatch("powermean: There must be one column weight " *
                                "per subcommunity"))
    @views map(axes(values, 2)) do col
        return _isemptycolumn(colweights, col) ? _emptymean(FP, orders) :
               powermean(values[:, col], orders, weights[:, col])
    end
end

"""
    qD

Calculates Hill / naive-similarity diversity of order(s) *qs* of a
population with given relative proportions.

# Arguments:

- `proportions`: relative proportions of different types in population

- `qs`: single number or vector of orders of diversity measurement

# Returns:

- Diversity of order qs (single number or vector of diversities)

"""
function qD end
function qD(asm::AbstractAssemblage, qs)
    hassimilarity(asm) &&
        error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    countsubcommunities(asm) == 1 ||
        throw(DimensionMismatch("Can only calculate diversity of a single community"))

    return powermean(getabundance(asm), qs .- 1, getabundance(asm))[1] .^ -1
end

function qD(proportions::AbstractVector{<:Real}, qs)
    return qD(Metacommunity(proportions), qs)
end

"""
    qDZ

Calculates Leinster-Cobbold / similarity-sensitive diversity of >= 1
order(s) *qs* of a population with given relative *proportions*, and
similarity matrix *Z*.

# Arguments:

- `proportions`: relative proportions of different types in a population

- `qs`: single number or vector of orders of diversity measurement

- `Z`: similarity matrix

# Returns:

- Diversity of order qs (single number or vector of diversities)

    """
function qDZ end

function qDZ(asm::AbstractAssemblage, qs)
    countsubcommunities(asm) == 1 ||
        throw(DimensionMismatch("Can only calculate diversity of a single community"))

    return powermean(getordinariness!(asm), qs .- 1, getabundance(asm))[1] .^ -1
end

function qDZ(proportions::AbstractVector{<:Real}, qs,
             sim::AbstractTypes = UniqueTypes(size(proportions, 1)))
    return qDZ(Metacommunity(proportions, sim), qs)
end

function qDZ(proportions::AbstractVector{<:Real}, qs,
             Z::AbstractMatrix{<:AbstractFloat})
    return qDZ(Metacommunity(proportions, GeneralTypes(Z)), qs)
end
