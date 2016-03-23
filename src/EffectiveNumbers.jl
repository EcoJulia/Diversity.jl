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

function powermean{S <: AbstractFloat}(values::Matrix{S},
                                       orders::Vector{S},
                                       weights::Matrix{S} = ones(values))
    map(order -> map(col -> powermean(values[:,col], order, weights[:, col]),
                     1:size(values)[2]), orders)
end

## We need to handle lack of automatic promotion between ints and floats in Julia
function powermean{S <: Real,
                   T <: Real,
                   U <: Real}(values::Array{S},
                              order::T,
                              weights::Array{U} = ones(values))
    powermean(values * 1., order * 1., weights * 1.)
end

## We need to handle matrices as well as vectors
function powermean{S <: Real,
                   T <: Real,
                   U <: Real}(values::Array{S},
                              orders::Vector{T},
                              weights::Array{U} = ones(values))
    map(order -> powermean(values * 1., order * 1., weights * 1.), orders)
end

"""
### Calculates Hill / naive-similarity diversity

Calculates Hill number or naive diversity of order(s) *qs* of a
population with given relative proportions.

#### Arguments:
- `proportions`: relative proportions of different individuals /
               species in population or series of populations
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

## We need to handle matrices as well as vectors
qD{S <: Real}(proportions::Matrix{S}, qs) =
    mapslices(p -> qD(p, qs), proportions, 1)

"""
### Calculates Leinster-Cobbold / similarity-sensitive diversity

Calculates Leinster-Cobbold general diversity of >= 1 order(s) *qs* of
a population with given relative *proportions*, and similarity matrix
*Z*.

#### Arguments:
- `proportions`: relative proportions of different individuals /
               species in a population or series of populations
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

## We need to handle matrices as well as vectors
qDZ{S <: AbstractFloat}(proportions::Matrix{S}, qs,
                        Z::Matrix{S} = eye(size(proportions, 1))) =
    mapslices(p -> qDZ(p, qs, Z), proportions, 1)
