module TestMetacommunity
using Diversity
using DataFrames
using Base.Test

three = [0.3, 0.3, 0.4]
oc_count = Onecommunity()
@testset "Onecommunity" begin
    oc_2 = Onecommunity()
    @test Diversity.countsubcommunities(oc_count) == Diversity.countsubcommunities(oc_2)
end

sim = [1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 1.0 1.0]
ab3 = [1 2; 3 0; 0 4]'
abnorm = ab3 / sum(ab3)
ms = GeneralTypes(sim)
@testset "GeneralTypes" begin
    @test getsimilarity(ms) == sim
    @test_throws DomainError GeneralTypes(-sim)
    @test_throws MethodError GeneralTypes(ab3)
    @test_throws DimensionMismatch GeneralTypes(convert(Matrix{Float64}, ab3))
end

tax = Taxonomy(DataFrame(Species=["This", "That"]), Dict(:Species=>1.0))
@testset "Taxonomy / psmatch" begin
    @test_throws ErrorException Diversity.getsimilarity(tax)
    @test_throws ErrorException Diversity.getordinariness(tax, abnorm)
    @test_throws ErrorException Taxonomy(DataFrame(Species=["This", "That"]),
                                         Dict(:Species=>1.0, :Genus=>0.5))
end

sc = Subcommunities(size(ab3, 2))
meta = Metacommunity(three, ms, oc_count)
sp = Species(size(ab3, 1))
meta2 = Metacommunity(convert(Matrix{Float64}, ab3), sp, sc, true)
@testset "Metacommunity" begin
    @test gettypes(meta) == ms
    @test getpartition(meta) == oc_count
    @test isnull(meta.ordinariness)
    @test getordinariness!(meta) ≈ [0.3, 0.6, 1.0]
    @test !isnull(meta.ordinariness)
    @test_throws DimensionMismatch Metacommunity(ab3, ms, sc)
    @test getsimilarity(meta2) ≈ eye(size(ab3, 1))
end

end

