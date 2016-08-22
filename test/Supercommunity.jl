module TestSupercommunity
using Diversity
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

oc_count = Onecommunity([3, 3, 4])
@testset "Onecommunity" begin
    # counts (int) are normalised, relative abundances (float) should sum to 1.0
    oc_rel = Onecommunity([0.3, 0.3, 0.4])
    @test oc_count.abundances ≈ oc_rel.abundances
    # force normalisation
    oc_manual = Onecommunity([3.0, 3.0, 4.0], true)
    @test oc_count.abundances ≈ oc_manual.abundances
    # Not normalised
    @test_throws ErrorException Onecommunity([3.0, 3.0, 4.0])
end

sim = [1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 1.0 1.0]
ab3 = [1 2; 3 0; 0 4]'
ms = MatrixSimilarity(sim)
@testset "MatrixSimilarity" begin
    @test getsimilarity(oc_count, ms) ≈ sim
    @test_throws DomainError MatrixSimilarity(-sim)
    @test_throws MethodError MatrixSimilarity(ab3)
    @test_throws DimensionMismatch MatrixSimilarity(convert(Matrix{Float64},
                                                            ab3))
end

tax = Taxonomy(Dict{AbstractString,
               Tuple{Float64, Dict{AbstractString,
               AbstractString}}}())
@testset "Taxonomy / psmatch" begin
    @test_throws ErrorException Diversity.psmatch(oc_count, tax)
    @test_throws ErrorException Diversity.getsimilarity(oc_count, tax)
    @test_throws ErrorException Diversity.getordinariness(oc_count, tax)
    @test Diversity.psmatch(oc_count, Unique()) == true
end

sc = Subcommunities(ab3)
sup = Supercommunity(oc_count, ms)
sp = Species()
sup2 = Supercommunity(sc, sp)
@testset "Supercommunity" begin
    @test sup.similarity == ms
    @test sup.partition == oc_count
    @test isnull(sup.ordinariness)
    @test getordinariness!(sup) ≈ [0.3, 0.6, 1.0]
    @test !isnull(sup.ordinariness)
    @test_throws DimensionMismatch Supercommunity(sc, ms)
    @test getsimilarity(sup2) ≈ eye(size(ab3, 1))
end

end

