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
function powermean{S <: AbstractFloat}(values::Vector{S},
                                       order::S = 1.,
                                       weights::Vector{S} = ones(values))
    length(values) == length(weights) ||
    throw(DimensionMismatch("powermean: Weight and value vectors must be the same length"))
    
    # Normalise weights to sum to 1 (as per Rényi)
    proportions = weights / sum(weights)

    # Check whether all proportions are NaN - happens in normalisation when all
    # weights are zero in group. In that case we want to propagate the NaN
    if (all(isnan(proportions)))
        return(NaN)
    end
    
    # Extract values with non-zero weights
    present = filter(x -> !isapprox(x[1], 0.), zip(proportions, values))
    if (isinf(order))
      if (order > 0.) # +Inf -> Maximum
        reduce((a, b) -> a[2] > b[2] ? a : b, present)[2]
      else # -Inf -> Minimum
        reduce((a, b) -> a[2] < b[2] ? a : b, present)[2]
      end
    else
      if (isapprox(order, 0))
        mapreduce(pair -> pair[2] ^ pair[1], *, present)
      else
        mapreduce(pair -> pair[1] * pair[2] ^ order, +,
                  present) ^ (1. / order)
        end
    end
end

# This is the next most simple case - matrices with subcommunities, and an order
function powermean{S <: AbstractFloat}(values::Matrix{S},
                                       order::S,
                                       weights::Matrix{S} = ones(values))
    size(values) == size(weights) ||
    throw(DimensionMismatch("powermean: Weight and value matrixes must be the same size"))

    map(col -> powermean(values[:,col], order, weights[:, col]), 1:size(values)[2])
end

# This is the next most common case - matrices with subcommunities, and a vector of orders
function powermean{S <: AbstractFloat}(values::Matrix{S},
                                       orders::Vector{S},
                                       weights::Matrix{S} = ones(values))
    map(order -> powermean(values, order, weights), orders)
end

## We need to handle lack of automatic promotion between ints and floats in Julia
function powermean{S <: Real,
                   T <: Real,
                   U <: Real}(values::Array{S},
                              order::T,
                              weights::Array{U} = ones(values))
    (length(size(values)) <= 2 && size(values) == size(weights)) ||
    throw(DimensionMismatch("powermean: Value and weight matrices must match and be at most 2-dimensional"))

    # Must be at least Float32
    FType = promote_type(S, T, U, Float32)
    powermean(convert(Array{FType}, values),
              convert(FType, order),
              convert(Array{FType}, weights))
end

## We need to handle lack of automatic promotion between ints and floats in Julia
function powermean{S <: Real,
                   T <: Real,
                   U <: Real}(values::Array{S},
                              orders::Vector{T},
                              weights::Array{U} = ones(values))
    (length(size(values)) <= 2) ||
    throw(DimensionMismatch("powermean: Value matrix must be at most 2-dimensional"))

    # Must be at least Float32
    FType = promote_type(S, T, U, Float32)
    map(order -> powermean(convert(Array{FType}, values), order, convert(Array{FType}, weights)),
        convert(Array{FType}, orders))
end

"""
### Calculates Hill / naive-similarity diversity

Calculates Hill number or naive diversity of order(s) *qs* of a
population with given relative proportions.

#### Arguments:
- `proportions`: relative proportions of different types in population
- `qs`: single number or vector of orders of diversity measurement

#### Returns:
- Diversity of order qs (single number or vector of diversities)"""
function qD{S <: Real}(proportions::Vector{S}, qs)
    if !isapprox(sum(proportions), 1.)
        warn("qD: Population proportions don't sum to 1, fixing...")
        proportions /= sum(proportions)
    end
    powermean(proportions, qs - 1., proportions) .^ -1
end

"""
### Calculates Leinster-Cobbold / similarity-sensitive diversity

Calculates Leinster-Cobbold general diversity of >= 1 order(s) *qs* of
a population with given relative *proportions*, and similarity matrix
*Z*.

#### Arguments:
- `proportions`: relative proportions of different types in a population
- `qs`: single number or vector of orders of diversity measurement
- `Z`: similarity matrix

#### Returns:
- Diversity of order qs (single number or vector of diversities)"""
function qDZ{S <: AbstractFloat}(proportions::Vector{S}, qs,
                                 Z::Matrix{S} = eye(length(proportions)))
    if !isapprox(sum(proportions), 1.)
        warn("qDZ: Population proportions don't sum to 1, fixing...")
        proportions /= sum(proportions)
    end

    l = length(proportions)
    size(Z) == (l, l) ||
    error("qDZ: Similarity matrix size does not match species number")
    powermean(Z * proportions, qs - 1., proportions) .^ -1
end
