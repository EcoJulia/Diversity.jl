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

end
