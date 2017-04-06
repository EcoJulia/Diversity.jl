using Diversity.ShortNames

function normsubalpha(meta::AbstractMetacommunity, qs)
    subdiv(ᾱ(meta), qs)
end
function rawsubalpha(meta::AbstractMetacommunity, qs)
    subdiv(α(meta), qs)
end
function normsubbeta(meta::AbstractMetacommunity, qs)
    subdiv(β̄(meta), qs)
end
function rawsubbeta(meta::AbstractMetacommunity, qs)
    subdiv(β(meta), qs)
end
function normsubrho(meta::AbstractMetacommunity, qs)
    subdiv(ρ̄(meta), qs)
end
function rawsubrho(meta::AbstractMetacommunity, qs)
    subdiv(ρ(meta), qs)
end
function subgamma(meta::AbstractMetacommunity, qs)
    subdiv(Γ(meta), qs)
end
function normmetaalpha(meta::AbstractMetacommunity, qs)
    metadiv(ᾱ(meta), qs)
end
function rawmetaalpha(meta::AbstractMetacommunity, qs)
    metadiv(α(meta), qs)
end
function normmetabeta(meta::AbstractMetacommunity, qs)
    metadiv(β̄(meta), qs)
end
function rawmetabeta(meta::AbstractMetacommunity, qs)
    metadiv(β(meta), qs)
end
function normmetarho(meta::AbstractMetacommunity, qs)
    metadiv(ρ̄(meta), qs)
end
function rawmetarho(meta::AbstractMetacommunity, qs)
    metadiv(ρ(meta), qs)
end
function metagamma(meta::AbstractMetacommunity, qs)
    metadiv(Γ(meta), qs)
end

"""
### Calculates subcommunity and metacommunity diversities

Calculates any diversity of a Metacommunity for a series of orders,
repesented as one or a vector of qs.

#### Arguments:
- `dls`: an iterable collection of DiversityLevels
- `dms`: an iterable collection of DiversityMeasures
- `meta`: a Metacommunity
- `qs`: single number or vector of values of parameter q

#### Returns:

A vector containing all of the diversity levels of all of the requested diversities.
"""
function diversity(dls, dms, meta::AbstractMetacommunity, qs)
    return mapreduce(measure -> mapreduce(dl -> dl(measure, qs), append!, dls),
                     append!, map(dm -> dm(meta), dms))
end
