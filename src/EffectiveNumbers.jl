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
function powermean(values::V, order::R = 1, weights::V = ones(values)) where
    {R <: Real, FP <: AbstractFloat, V <: AbstractVector{FP}}
    length(values) == length(weights) ||
    throw(DimensionMismatch("powermean: Weight and value vectors must be the same length"))

    # Normalise weights to sum to 1 (as per Rényi)
    proportions = weights / sum(weights)

    # Check whether all proportions are NaN - happens in normalisation when all
    # weights are zero in group. In that case we want to propagate the NaN
    if all(isnan.(proportions))
        return(convert(FP, NaN))
    end

    # Extract values with non-zero weights
    present = [(proportions[i], values[i])
               for i=1:length(proportions) if !(proportions[i] ≈ zero(eltype(proportions)))]
    if isinf(order)
        if order > 0 # +Inf -> Maximum
            reduce((a, b) -> a[2] > b[2] ? a : b, present)[2]
        else # -Inf -> Minimum
            reduce((a, b) -> a[2] < b[2] ? a : b, present)[2]
        end
    else
        if order ≈ zero(order)
            mapreduce(pair -> pair[2] ^ pair[1], *, present)
        else
            mapreduce(pair -> pair[1] * pair[2] ^ order, +, present) ^
                (one(order) / order)
        end
    end
end

# This is the next most common case - a vector of orders
function powermean(values::V,
                   orders::VR,
                   weights::V = ones(values)) where
    {R <: Real, VR <: AbstractVector{R},
     FP <: AbstractFloat, V <: AbstractVector{FP}}
    map(order -> powermean(values, order, weights), orders)
end

# This is the next most simple case - matrices with subcommunities, and an order or orders
function powermean(values::M, orders, weights::M = ones(values)) where
    {FP <: AbstractFloat, M <: AbstractMatrix{FP}}
    size(values) == size(weights) ||
        throw(DimensionMismatch("powermean: Weight and value matrixes " *
                                "must be the same size"))
    
    map(col -> powermean(view(values, :, col), orders,
                         view(weights, :, col)), 1:size(values, 2))
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

function qD(meta::M, qs) where M <: AbstractMetacommunity
    countsubcommunities(meta) == 1 ||
    throw(DimensionMismatch("Can only calculate diversity of a single community"))

    powermean(getabundance(meta), qs - 1, getabundance(meta))[1] .^ -1
end

function qD(proportions::V, qs) where {FP <: AbstractFloat,
                                       V <: AbstractVector{FP}}
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

function qDZ(meta::M, qs) where M <: AbstractMetacommunity
    countsubcommunities(meta) == 1 ||
    throw(DimensionMismatch("Can only calculate diversity of a single community"))

    powermean(getordinariness!(meta), qs - 1, getabundance(meta))[1] .^ -1
end

function qDZ(proportions::V, qs,
             sim::Sim = UniqueTypes(size(proportions, 1))) where
    {FP <: AbstractFloat, V <: AbstractVector{FP}, Sim <: AbstractTypes}
    qDZ(Metacommunity(proportions, sim), qs)
end

function qDZ(proportions::V, qs, Z::M) where
    {FP <: AbstractFloat, V <: AbstractVector{FP}, M <: AbstractMatrix{FP}}
    qDZ(Metacommunity(proportions, GeneralTypes(Z)), qs)
end
