## contributions - Calculate diversity contributions from sub-communities
## Calculates 
##
## Arguments:
## - proportions - population proportions
## - measure - diversity measure to use
## - perindividual - do we measure per individual in population (true)
##                   or per sub-community (false)
## - Z - similarity matrix
##
## Returns:
## - contributions
function contributions{S <: FloatingPoint}(proportions::Matrix{S}, qs,
                                           measure::Function,
                                           perindividual::Bool = true,
                                           Z::Matrix{S} = eye(size(proportions)[1]))
    diversities = measure(proportions, qs, Z)
    indices = 1:length(qs)
    weights = mapslices(sum, proportions, 1)
    results = 0. .* diversities
    for (i in indices)
        power = 1. - qs[i]
        if (isinf(power))
            if (power > 0) # +Inf -> Maximum, so 1 for top otherwise 0
                mx = 0. * diversities[i, :]
                mx[indmax(diversities[i, :])] = 1.
                results[i, :] = mx
            else # -Inf -> Minimum, so 1 for bottom otherwise 0
                mn = 0. * diversities[i, :]
                mn[indmin(diversities[i, :])] = 1.
                results[i, :] = mn
            end
        else
            ## We want to use scaled values so that ultimately we can
            ## get everything relative to total and find out
            ## contribution per individual or community rather than
            ## any kind of actual value
            if (isapprox(power, 0))
                results[i, :] = weights .* log(diversities[i, :])
                ## IS THIS THE RIGHT WAY TO NORMALISE UNDER LOG?
                results[i, :] -= weights .* sum(results[i, :])
            else
                results[i, :] = weights .* (diversities[i, :] .^ power)
                results[i, :] /= sum(results[i, :])
            end
        end
        
    end
    perindividual ? (results ./ weights) : results
end
