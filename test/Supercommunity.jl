module TestSupercommunity
using Diversity
using Base.Test
using Compat

oc = Onecommunity([3.0, 3.0, 4.0], true)
@test_approx_eq Onecommunity([0.3, 0.3, 0.4]).abundances oc.abundances

sim = [1.0 0 0; 1.0 1.0 0.0; 1.0 1.0 1.0]
ms = MatrixSimilarity(sim)
@test_approx_eq getSimilarityMatrix(oc, ms) sim
@test_throws DomainError MatrixSimilarity(-sim)

sup = Supercommunity(oc, ms)
@test sup.similarity == ms
@test sup.partition == oc
@test isnull(sup.ordinariness)

@test_approx_eq getOrdinariness!(sup) [0.3, 0.6, 1.0]
@test !isnull(sup.ordinariness)

tax = Taxonomy(Dict{AbstractString, @compat(Tuple{Float64, Dict{AbstractString, AbstractString}})}())
@test_throws ErrorException Diversity.match!(oc, tax)

ab3 = [1.0 2.0; 3.0 0.0; 0.0 4.0]'
@test_throws DimensionMismatch MatrixSimilarity(ab3)
sc = Subcommunities(ab3, true)
@test_throws DimensionMismatch Supercommunity(sc, ms)

sp = Species()
sup2 = Supercommunity(sc, sp)
@test_approx_eq getSimilarityMatrix(sup2) eye(size(ab3, 1))

end

