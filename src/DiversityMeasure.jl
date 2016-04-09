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
