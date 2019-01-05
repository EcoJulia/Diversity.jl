using EcoBase: AbstractAssemblage
using Diversity.ShortNames

function norm_sub_alpha(meta::AbstractAssemblage, qs)
    subdiv(ᾱ(meta), qs)
end
function raw_sub_alpha(meta::AbstractAssemblage, qs)
    subdiv(α(meta), qs)
end
function norm_sub_beta(meta::AbstractAssemblage, qs)
    subdiv(β̄(meta), qs)
end
function raw_sub_beta(meta::AbstractAssemblage, qs)
    subdiv(β(meta), qs)
end
function norm_sub_rho(meta::AbstractAssemblage, qs)
    subdiv(ρ̄(meta), qs)
end
function raw_sub_rho(meta::AbstractAssemblage, qs)
    subdiv(ρ(meta), qs)
end
function sub_gamma(meta::AbstractAssemblage, qs)
    subdiv(Γ(meta), qs)
end
function norm_meta_alpha(meta::AbstractAssemblage, qs)
    metadiv(ᾱ(meta), qs)
end
function raw_meta_alpha(meta::AbstractAssemblage, qs)
    metadiv(α(meta), qs)
end
function norm_meta_beta(meta::AbstractAssemblage, qs)
    metadiv(β̄(meta), qs)
end
function raw_meta_beta(meta::AbstractAssemblage, qs)
    metadiv(β(meta), qs)
end
function norm_meta_rho(meta::AbstractAssemblage, qs)
    metadiv(ρ̄(meta), qs)
end
function raw_meta_rho(meta::AbstractAssemblage, qs)
    metadiv(ρ(meta), qs)
end
function meta_gamma(meta::AbstractAssemblage, qs)
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
function diversity(dls, dms, meta::AbstractAssemblage, qs)
    return mapreduce(measure -> mapreduce(dl -> dl(measure, qs), append!, dls),
                     append!, map(dm -> dm(meta), dms))
end
