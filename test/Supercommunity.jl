module TestSupercommunity
using Diversity
using Base.Test

oc = Onecommunity([3.0, 3.0, 4.0], true)
#@testset "OneCommunity" begin
@test_approx_eq Onecommunity([0.3, 0.3, 0.4]).abundances oc.abundances
#end

sim = [1.0 0 0; 1.0 1.0 0.0; 1.0 1.0 1.0]
ab3 = [1.0 2.0; 3.0 0.0; 0.0 4.0]'
ms = MatrixSimilarity(sim)
#@testset "MatrixSimilarity" begin
@test_approx_eq getsimilarity(oc, ms) sim
@test_throws DomainError MatrixSimilarity(-sim)
@test_throws DimensionMismatch MatrixSimilarity(ab3)
#end

tax = Taxonomy(Dict{AbstractString,
               Tuple{Float64, Dict{AbstractString,
               AbstractString}}}())
#@testset "Taxonomy" begin
@test_throws ErrorException Diversity.psmatch(oc, tax)
#end

sc = Subcommunities(ab3, true)
sup = Supercommunity(oc, ms)
sp = Species()
sup2 = Supercommunity(sc, sp)
#@testset "Supercommunity" begin
@test sup.similarity == ms
@test sup.partition == oc
@test isnull(sup.ordinariness)
@test_approx_eq getordinariness!(sup) [0.3, 0.6, 1.0]
@test !isnull(sup.ordinariness)
@test_throws DimensionMismatch Supercommunity(sc, ms)
@test_approx_eq getsimilarity(sup2) eye(size(ab3, 1))
#end

end

