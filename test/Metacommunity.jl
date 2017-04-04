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
@testset "Taxonomy" begin
    @test Diversity.counttypes(tax) == 2
    @test Diversity.hasnames(tax)
    @test Diversity.getnames(tax) == ["This", "That"]
    @test_throws ErrorException Diversity.getsimilarity(tax)
    @test_throws ErrorException Diversity.getordinariness(tax, abnorm)
    @test_throws ErrorException Taxonomy(DataFrame(Species=["This", "That"]),
                                         Dict(:Species=>1.0, :Genus=>0.5))
end

@testset "Type names" begin
    @test !Diversity.hasnames(Species(3))
    @test Diversity.hasnames(Species(["My species"]))
    @test Diversity.getnames(GeneralTypes(eye(1), ["My species"])) == ["My species"]
    @test Diversity.getnames(UniqueTypes(["One", "Two"])) == ["One", "Two"]
end

sc = Subcommunities(size(ab3, 2))
meta = Metacommunity(three, ms, oc_count)
sp = Species(size(ab3, 1))
abf = ab3 ./ sum(ab3)
meta2 = Metacommunity(abf, sp, sc)
@testset "Metacommunity" begin
    @test Diversity.gettypes(meta) == ms
    @test Diversity.getpartition(meta) == oc_count
    @test isnull(meta.ordinariness)
    @test getordinariness!(meta) ≈ [0.3, 0.6, 1.0]
    @test !isnull(meta.ordinariness)
    @test_throws TypeError Metacommunity(ab3, ms, sc)
    @test_throws ErrorException Metacommunity(-abf, ms, sc)
    @test Diversity.getsimilarity(Diversity.gettypes(meta2)) ≈ eye(size(ab3, 1))
end

end

