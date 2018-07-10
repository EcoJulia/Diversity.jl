using Diversity

"""
    hillnumber(proportions, qs)

Calculate the Hill number (or naive diversity) of order q of
population(s) with given relative proportions

# Arguments:

- `proportions`: relative proportions of different individuals / species
                 in population (vector, or matrix where columns are
                 individual populations)

- `qs`: single number or vector of orders of diversity measurement

# Returns:

- Diversity of order qs (single number or vector of diversities)
"""
function hillnumber(proportions, qs)
    hill = subdiv(NormalisedAlpha(Metacommunity(proportions)), qs)
    hill[:measure] = "HillNumber"
    delete!(hill, :div_type)
    return hill
end
