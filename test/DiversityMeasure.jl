module TestDiversityMeasure
using Diversity
using Base.Test
using Compat

oc = Onecommunity([3.0, 3.0, 4.0], true)
sim = [1.0 0 0; 1.0 1.0 0.0; 1.0 1.0 1.0]
ms = MatrixSimilarity(sim)
sup = Supercommunity(oc, ms)

ab3 = [1.0 2.0; 3.0 0.0; 0.0 4.0]'
sc = Subcommunities(ab3, true)
sp = Species()
sup2 = Supercommunity(sc, sp)

nr = NormalisedRho(sup)
@test getName(nr) == "ρ̄"
@test getASCIIName(nr) == "rho bar"
@test getFullName(nr) == "representativeness"

sup = Supercommunity(oc)
nab = NormalisedAlpha(sup2)

@test getName(nab) == "ᾱ"
@test getASCIIName(nab) == "alpha bar"
@test getFullName(nab) == "normalised alpha diversity"

scg = subcommunityDiversity(Gamma(sup2))
#scg([0.0, 1.0, 2.0])

#individualDiversity(nab, 0.0)
#subcommunityDiversity(nab, 0.0)

end
