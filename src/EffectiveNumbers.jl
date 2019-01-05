"""
### Calculates the weighted powermean of a series of numbers

Calculates *order*th power mean of *values*, weighted by
*weights*. By default, *weights* are equal and *order*
is 1, so this is just the arithmetic mean.

#### Arguments:
- `values`: values for which to calculate mean
- `order`: order of power mean
- `weights`: weights of elements, normalised to 1 inside function

#### Returns:
- weighted power mean(s)
"""
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
                    s += weights[i] * values[i] ^ order
                end
            end
            return (s / sum(weights)) ^ (one(FP) / order)
        end
    end
end

# This is the next most common case - a vector of orders
function powermean(values::V1,
                   orders::VR,
                   weights::V2 = fill!(similar(values), 1.0)) where
    {R <: Real, VR <: AbstractVector{R},
     FP <: AbstractFloat, V1 <: AbstractVector{FP}, V2 <: AbstractVector{FP}}
    return map(order -> powermean(values, order, weights), orders)
end

# This is the next most simple case - matrices with subcommunities, and an order or orders
function powermean(values::M1, orders, weights::M2 = fill!(similar(values), 1)) where
    {FP <: AbstractFloat, M1 <: AbstractMatrix{FP}, M2 <: AbstractMatrix{FP}}
    size(values) == size(weights) ||
        throw(DimensionMismatch("powermean: Weight and value matrixes " *
                                "must be the same size"))
    @views map(col -> powermean(values[:, col], orders,
                                weights[:, col]), 1:size(values, 2))
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
    hassimilarity(asm) && error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    countsubcommunities(asm) == 1 ||
    throw(DimensionMismatch("Can only calculate diversity of a single community"))

    powermean(getabundance(asm), qs .- 1, getabundance(asm))[1] .^ -1
end

function qD(proportions::AbstractVector{<: Real}, qs)
    qD(Metacommunity(proportions), qs)
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

    powermean(getordinariness!(asm), qs .- 1, getabundance(asm))[1] .^ -1
end

function qDZ(proportions::AbstractVector{<: Real}, qs,
             sim::AbstractTypes = UniqueTypes(size(proportions, 1)))
    qDZ(Metacommunity(proportions, sim), qs)
end

function qDZ(proportions::AbstractVector{<: Real}, qs, Z::AbstractMatrix{<: AbstractFloat})
    qDZ(Metacommunity(proportions, GeneralTypes(Z)), qs)
end
