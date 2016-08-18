module TestDiversityMeasure
using Diversity
using Base.Test
using Compat

oc = Onecommunity([3, 3, 4])
sim = [1.0 0 0; 1.0 1.0 0.0; 1.0 1.0 1.0]
ms = MatrixSimilarity(sim)
sup = Supercommunity(oc, ms)
sup3 = Supercommunity(oc)

ab3 = [1 2; 3 0; 0 4]'
sc = Subcommunities(ab3)
sp = Species()
sup2 = Supercommunity(sc, sp)
scg = subcommunityDiversity(Gamma(sup2))
@test_approx_eq scg(1) scg(1.0)

nr = NormalisedRho(sup)
@test getName(nr) == "ρ̄"
@test getASCIIName(nr) == "rho bar"
@test getFullName(nr) == "representativeness"

nab = NormalisedAlpha(sup2)

@test getName(nab) == "ᾱ"
@test getASCIIName(nab) == "alpha bar"
@test getFullName(nab) == "normalised alpha diversity"

individualDiversity(nab, 0.0)
subcommunityDiversity(nab, 0.0)

end
