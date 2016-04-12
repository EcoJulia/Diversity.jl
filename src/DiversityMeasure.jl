"""
### Generates the function to calculate individual diversities

Generates the function to calculate individual diversities for a
series of orders, represented as a vector of qs.

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- Function which takes a single number or vector of values of
  parameter q, and returns the individual diversities for those
  values.
"""
:individualDiversity

"""
### Generates the function to calculate subcommunity diversity

Generates the function to calculate subcommunity diversity for a
series of orders, represented as a vector of qs.

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- Function which takes a single number or vector of values of
  parameter q, and returns the subcommunity diversities for those values.
"""
:subcommunityDiversity

"""
### Generates the function to calculate supercommunity diversity

Generates the function to calculate supercommunity diversity for a
series of orders, represented as a vector of qs.

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- Function which takes a single number or vector of values of
  parameter q, and returns the supercommunity diversities for those
  values.
"""
:supercommunityDiversity

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

"""
### Return the name of the DiversityMeasure

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- String containing true unicode name of DiversityMeasure
"""
function getName(dm::DiversityMeasure)
    replace(string(typeof(dm)), "Diversity.", "")
end

"""
### Return the ASCII name of the DiversityMeasure

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- String containing simple ASCII name of DiversityMeasure
"""
function getASCIIName
end

"""
### Return the full name of the DiversityMeasure

#### Arguments:

- `dm`: DiversityMeasure

#### Returns:

- String containing full descriptive name of DiversityMeasure
"""
function getFullName
end

call(dl::DiversityLevel, dm::DiversityMeasure) = getPartitionFunction(dm, dl)
call(dl::DiversityLevel, dm::DiversityMeasure, others...) = getPartitionFunction(dm, dl)(others...)

"""
### Supertype of all power mean-based diversity measures

This abstract DiversityMeasure subtype is the supertype of all
diversity measures which are straight power means. PowerMeanMeasure
subtypes allow you to calculate and cache any kind of diversity of a
supercommunity.
"""
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

"""
### Supertype of all relative entropy-based diversity measures

This abstract DiversityMeasure subtype is the supertype of all
diversity measures which are straight power means.
RelativeEntropyMeasure subtypes allow you to calculate and cache any
kind of diversity of a supercommunity.
"""
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

"""
### Raw alpha diversity type (α)

Calculates raw alpha diversity (α) of all of the individuals in a
supercommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

#### Constructor arguments:

- `sup`: a Supercommunity
"""
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

"""
### Normalised alpha diversity type (ᾱ)

Calculates normalised alpha diversity (ᾱ) of all of the individuals in
a supercommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

#### Constructor arguments:

- `sup`: a Supercommunity
"""
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

"""
### Distinctiveness (β, raw beta diversity) type

Calculates distinctiveness (β, raw beta diversity) of all of the individuals in a
supercommunity, and caches them for subsequent analysis. This is a
subtype of RelativeEntropyMeasure, meaning that subcommunity and type
composite diversity measures are relative entropies, and their
composite types are powermeans of those measures.

#### Constructor arguments:

- `sup`: a Supercommunity
"""
type β <: RelativeEntropyMeasure
    abundances::Array
    diversities::Array
    function ρ(sup::AbstractSupercommunity)
        new(getAbundances(sup),
            getSuperOrdinariness!(sup) ./ getOrdinariness!(sup))
    end
end

typealias RawBeta β
typealias Distinctiveness β

getASCIIName(::β) = "beta"
getFullName(::β) = "distinctiveness"

"""
### Normalised beta diversity type (β̄)

Calculates normalised beta diversity (β̄) of all of the individuals in
a supercommunity, and caches them for subsequent analysis. This is a
subtype of RelativeEntropyMeasure, meaning that subcommunity and type
composite diversity measures are relative entropies, and their
composite types are powermeans of those measures.

#### Constructor arguments:

- `sup`: a Supercommunity
"""
type β̄ <: RelativeEntropyMeasure
    abundances::Array
    diversities::Array
    function β̄(sup::AbstractSupercommunity)
        new(getAbundances(sup),
            getOrdinariness!(sup) ./ (getSuperOrdinariness!(sup) * getWeights(sup)))
    end
end

typealias NormalisedBeta β̄

getASCIIName(::β̄) = "beta bar"
getFullName(::β̄) = "effective number of subcommunities"

"""
### Redundancy (ρ, raw beta diversity) type

Calculates redundancy (ρ, raw beta diversity) of all of the
individuals in a supercommunity, and caches them for subsequent
analysis. This is a subtype of PowerMeanMeasure, meaning that all
composite diversity measures are simple powermeans of the individual
measures.

#### Constructor arguments:

- `sup`: a Supercommunity
"""
type ρ <: PowerMeanMeasure
    abundances::Array
    diversities::Array
    function ρ(sup::AbstractSupercommunity)
        new(getAbundances(sup),
            getSuperOrdinariness!(sup) ./ getOrdinariness!(sup))
    end
end

typealias RawRho ρ
typealias Redundancy ρ

getASCIIName(::ρ) = "rho"
getFullName(::ρ) = "redundancy"

"""
### Representativeness (ρ̄, normalised beta diversity) type

Calculates redundancy (ρ̄, normalised beta diversity) of all of the
individuals in a supercommunity, and caches them for subsequent
analysis. This is a subtype of PowerMeanMeasure, meaning that all
composite diversity measures are simple powermeans of the individual
measures.

#### Constructor arguments:

- `sup`: a Supercommunity
"""
type ρ̄ <: PowerMeanMeasure
    abundances::Array
    diversities::Array
    function ρ̄(sup::AbstractSupercommunity)
        new(getAbundances(sup),
            getSuperOrdinariness!(sup) * getWeights(sup) ./ getOrdinariness!(sup))
    end
end

typealias NormalisedRho ρ̄
typealias Representativeness ρ̄

getASCIIName(::ρ̄) = "rho bar"
getFullName(::ρ̄) = "representativeness"

"""
### Gamma diversity type (γ)

Calculates gamma diversity (γ) of all of the individuals in a
supercommunity, and caches them for subsequent analysis. This is a
subtype of PowerMeanMeasure, meaning that all composite diversity
measures are simple powermeans of the individual measures.

#### Constructor arguments:

- `sup`: a Supercommunity
"""
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
