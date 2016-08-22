"""
### Calculates subcommunity and supercommunity diversities

Calculates any diversity of a Supercommunity for a series of orders,
repesented as one or a vector of qs.

#### Arguments:
- `dl`: a DiversityLevel
- `dm`: a DiversityMeasure
- `sup`: a Supercommunity
- `qs`: single number or vector of values of parameter q

#### Returns:

The requested diversities.
"""
function diversity{Sup <: AbstractSupercommunity}(dl::DiversityLevel,
                                                  dm, sup::Sup, qs)
    dl(dm(sup), qs)
end

"""
### Calculates subcommunity and supercommunity diversities

Calculates any diversity of a Supercommunity for a series of orders,
repesented as one or a vector of qs.

#### Arguments:
- `dls`: a Set of DiversityLevels
- `dms`: a Set of DiversityMeasures
- `sup`: a Supercommunity
- `qs`: single number or vector of values of parameter q

#### Returns:

A vector containing all of the diversity levels of all of the requested diversities.
"""
function diversity{Sup <: AbstractSupercommunity}(dls::Set{DiversityLevel},
                                                  dms::Set, sup::Sup, qs)
    ret = Vector()
    for dm in dms
        dmv = dm(sup)
        for dl in dls
            push!(ret, dl(dmv, qs))
        end
    end
    ret
end
