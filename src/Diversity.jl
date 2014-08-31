module Diversity

export powermean

## powermean - Calculate order-th power mean of values, weighted by weights
## By default, weights are equal and order is 1, so this is just the arithmetic mean
##
## Arguments:
## - values - values for which to calculate mean
## - order - order of power mean
## - weights - weights of elements, normalised to 1 inside function
##
## Returns:
## - weighted power mean
function powermean(values::Array,
                   order = 1,
                   weights::Array = ones(FloatingPoint, size(values)))
    ## Normalise weights to sum to 1 (as per Rényi)
    proportions = weights / sum(weights)
    power = convert(FloatingPoint, order)
    present = filter(x -> !isapprox(x[1], 0), zip(proportions, values))
    if (isinf(power))
        if (power > 0) # Maximum
            println(reduce((a, b) -> a[2] > b[2] ? a : b, present)[2])
        else # Minimum
            println(reduce((a, b) -> a[2] < b[2] ? a : b, present)[2])
        end
    else
        if (isapprox(power, 0))
            mapreduce((pair) -> pair[2] ^ pair[1], *, present)
        else
            mapreduce(pair -> pair[1] * pair[2] ^ power, +,
                      present) ^ (1 / power)
        end 
    end
end

end # module
