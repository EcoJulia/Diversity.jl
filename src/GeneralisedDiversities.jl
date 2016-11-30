function normsubalpha(meta::AbstractMetacommunity, qs)
    subcommunityDiversity(ᾱ(meta), qs)
end
function rawsubalpha(meta::AbstractMetacommunity, qs)
    subcommunityDiversity(α(meta), qs)
end
function normsubbeta(meta::AbstractMetacommunity, qs)
    subcommunityDiversity(β̄(meta), qs)
end
function rawsubbeta(meta::AbstractMetacommunity, qs)
    subcommunityDiversity(β(meta), qs)
end
function normsubrho(meta::AbstractMetacommunity, qs)
    subcommunityDiversity(ρ̄(meta), qs)
end
function rawsubrho(meta::AbstractMetacommunity, qs)
    subcommunityDiversity(ρ(meta), qs)
end
function subgamma(meta::AbstractMetacommunity, qs)
    subcommunityDiversity(γ(meta), qs)
end
function normmetaalpha(meta::AbstractMetacommunity, qs)
    metacommunityDiversity(ᾱ(meta), qs)
end
function rawmetaalpha(meta::AbstractMetacommunity, qs)
    metacommunityDiversity(α(meta), qs)
end
function normmetabeta(meta::AbstractMetacommunity, qs)
    metacommunityDiversity(β̄(meta), qs)
end
function rawmetabeta(meta::AbstractMetacommunity, qs)
    metacommunityDiversity(β(meta), qs)
end
function normmetarho(meta::AbstractMetacommunity, qs)
    metacommunityDiversity(ρ̄(meta), qs)
end
function rawmetarho(meta::AbstractMetacommunity, qs)
    metacommunityDiversity(ρ(meta), qs)
end
function metagamma(meta::AbstractMetacommunity, qs)
    metacommunityDiversity(γ(meta), qs)
end

"""
### Calculates subcommunity and metacommunity diversities

Calculates any diversity of a Metacommunity for a series of orders,
repesented as one or a vector of qs.

#### Arguments:
- `dl`: a DiversityLevel
- `dm`: a DiversityMeasure
- `meta`: a Metacommunity
- `qs`: single number or vector of values of parameter q

#### Returns:

The requested diversities.
"""
function diversity{Meta <: AbstractMetacommunity}(dl::DiversityLevel,
                                                  dm, meta::Meta, qs)
    dl(dm(meta), qs)
end

"""
### Calculates subcommunity and metacommunity diversities

Calculates any diversity of a Metacommunity for a series of orders,
repesented as one or a vector of qs.

#### Arguments:
- `dls`: a Set of DiversityLevels
- `dms`: a Set of DiversityMeasures
- `meta`: a Metacommunity
- `qs`: single number or vector of values of parameter q

#### Returns:

A vector containing all of the diversity levels of all of the requested diversities.
"""
function diversity{Meta <: AbstractMetacommunity}(dls::Set{DiversityLevel},
                                                  dms::Set, meta::Meta, qs)
    ret = Vector()
    for dm in dms
        dmv = dm(meta)
        for dl in dls
            push!(ret, dl(dmv, qs))
        end
    end
    ret
end
