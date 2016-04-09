"""
### Enumeration of levels that can exist / be calculated for a supercommunity.
"""
@enum DiversityLevel individualDiversity subcommunityDiversity communityDiversity typeDiversity typeCollectionDiversity supercommunityDiversity metacommunityDiversity

"""
### DiversityMeasure supertype for all diversity measure types

This type is the abstract superclass of all diversity measure types.
DiversityMeasure subtypes allow you to calculate and cache any kind of
diversity of a supercommunity.
"""
abstract DiversityMeasure

function getName(div::DiversityMeasure)
    replace(string(typeof(div)), "Diversity.", "")
end

call(dl::DiversityLevel, dm::DiversityMeasure) = getPartitionFunction(dm, dl)

abstract PowerMeanMeasure <: DiversityMeasure
function getPartitionFunction(measure::PowerMeanMeasure,
                              level::DiversityLevel)
    if (level == individualDiversity)
        return function (qs)
            map(q -> measure.diversities, qs)
        end
    elseif (level == subcommunityDiversity)
        function (qs)
            divAllFn = getPartitionFunction(measure, individualDiversity)
            abundances = measure.abundances
            println(abundances)
            println(1.0 - qs)
            println(divAllFn(qs))
            map((order, divAll) -> {println(order);
                                    println(divAll);
                                    powermean(divAll, order, measure.abundances)},
                1.0 - qs, divAllFn(qs))
        end
    elseif (level == supercommunityDiversity)
        function (qs)
            divSubFn = getPartitionFunction(measure, subcommunityDiversity)
            map((order, divSub) -> powermean(divSub, order, measure.weights),
                1.0 - qs, divSubFn(qs))
        end
    else
        error("unrecognised request")
    end
end

abstract RelativeEntropyMeasure <: DiversityMeasure
function getPartitionFunction(measure::RelativeEntropyMeasure,
                              level::DiversityLevel)
    if (level == individualDiversity)
        return function (qs)
            map(q -> measure.diversities, qs)
        end
    elseif (level == subcommunityDiversity)
        function (qs)
            divAllFn = getPartitionFunction(measure, individualDiversity)
            abundances = measure.abundances
            println(abundances)
            println(1.0 - qs)
            println(divAllFn(qs))
            map((order, divAll) -> {println(order);
                                    println(divAll);
                                    powermean(divAll, order, measure.abundances)},
                qs - 1.0, divAllFn(qs))
        end
    elseif (level == supercommunityDiversity)
        function (qs)
            divSubFn = getPartitionFunction(measure, subcommunityDiversity)
            map((order, divSub) -> powermean(divSub, order, measure.weights),
                1.0 - qs, divSubFn(qs))
        end
    else
        error("unrecognised request")
    end
end
