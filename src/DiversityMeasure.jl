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

type α <: PowerMeanMeasure
    abundances::Array
    diversities::Array
    function α(sup::AbstractSupercommunity)
        new(getAbundances(sup),
            getOrdinariness!(sup) .^ -1)
    end
end

typealias RawAlpha α

getASCIIName(::α) = "alpha"
getFullName(::α) = "raw alpha diversity"

type ᾱ <: PowerMeanMeasure
    abundances::Array
    diversities::Array
    function ᾱ(sup::AbstractSupercommunity)
        new(getAbundances(sup),
            getWeights(sup) ./ getOrdinariness!(sup))
    end
end

typealias NormalisedAlpha ᾱ

getASCIIName(::ᾱ) = "alpha bar"
getFullName(::ᾱ) = "normalised alpha diversity"

type β <: RelativeEntropyMeasure
    abundances::Array
    diversities::Array
    function ρ(sup::AbstractSupercommunity)
        new(getAbundances(sup),
            getSuperOrdinariness!(sup) ./ getOrdinariness!(sup))
    end
end

typealias RawBeta β

getASCIIName(::β) = "beta"
getFullName(::β) = "distinctiveness"

type β̄ <: RelativeEntropyMeasure
    abundances::Array
    diversities::Array
    function ρ̄(sup::AbstractSupercommunity)
        new(getAbundances(sup),
            getOrdinariness!(sup) ./ (getSuperOrdinariness!(sup) * getWeights(sup)))
    end
end

typealias NormalisedBeta β̄

getASCIIName(::β̄) = "beta bar"
getFullName(::β̄) = "effective number of subcommunities"

type ρ <: PowerMeanMeasure
    abundances::Array
    diversities::Array
    function ρ(sup::AbstractSupercommunity)
        new(getAbundances(sup),
            getSuperOrdinariness!(sup) ./ getOrdinariness!(sup))
    end
end

typealias RawRho ρ

getASCIIName(::ρ) = "rho"
getFullName(::ρ) = "redundancy"

type ρ̄ <: PowerMeanMeasure
    abundances::Array
    diversities::Array
    function ρ̄(sup::AbstractSupercommunity)
        new(getAbundances(sup),
            getSuperOrdinariness!(sup) * getWeights(sup) ./ getOrdinariness!(sup))
    end
end

typealias NormalisedRho ρ̄

getASCIIName(::ρ̄) = "rho bar"
getFullName(::ρ̄) = "representativeness"

type γ <: PowerMeanMeasure
    abundances::Array
    diversities::Array
    function γ(sup::AbstractSupercommunity)
        new(getAbundances(sup),
            ones(1, length(sup)) ./ getSuperOrdinariness!(sup))
    end
end

typealias Gamma γ

getASCIIName(::γ) = "gamma"
getFullName(::γ) = "gamma diversity"

#type γ̄ <: PowerMeanMeasure
#    abundances::Array
#    diversities::Array
#    function γ̄(sup::AbstractSupercommunity)
#        new(getAbundances(sup),
#            sum(getWeights(sup)) * ones(1, length(sup)) ./ getSuperOrdinariness!(sup))
#    end
#end

#typealias NormalisedGammaDiversity γ̄

# getASCIIName(::γ̄) = "gamma bar"
# getFullName(::γ̄) = "normalised gamma diversity"
