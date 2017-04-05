module TestMetacommunity
using Diversity
using DataFrames
using Base.Test

three = [0.3, 0.3, 0.4]
oc_count = Onecommunity()
ab3 = [1 2; 3 0; 0 4]'
abnorm = ab3 / sum(ab3)
sc = Subcommunities(size(ab3, 2))
@testset "Communities" begin
    oc_2 = Onecommunity()
    @test countsubcommunities(oc_count) == countsubcommunities(oc_2)
    @test !hasnames(oc_2)
    @test_throws NullException getnames(sc)
end

sim = [1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 1.0 1.0]
ms = GeneralTypes(sim)
@testset "GeneralTypes" begin
    @test getsimilarity(ms) == sim
    @test_throws DomainError GeneralTypes(-sim)
    @test_throws MethodError GeneralTypes(ab3)
    @test_throws DimensionMismatch GeneralTypes(convert(Matrix{Float64}, ab3))
end

tax = Taxonomy(DataFrame(Species=["This", "That"]), Dict(:Species=>1.0))
@testset "Taxonomy" begin
    @test counttypes(tax) == 2
    @test hasnames(tax)
    @test getnames(tax) == ["This", "That"]
    @test_throws ErrorException getsimilarity(tax)
    @test_throws ErrorException Diversity.getordinariness(tax, abnorm)
    @test_throws ErrorException Taxonomy(DataFrame(Species=["This", "That"]),
                                         Dict(:Species=>1.0, :Genus=>0.5))
end

@testset "Type names" begin
    @test !hasnames(Species(3))
    @test hasnames(Species(["My species"]))
    @test getnames(GeneralTypes(eye(1), ["My species"])) == ["My species"]
    @test getnames(UniqueTypes(["One", "Two"])) == ["One", "Two"]
end

meta = Metacommunity(three, ms, oc_count)
sp = Species(size(ab3, 1))
abf = ab3 ./ sum(ab3)
meta2 = Metacommunity(abf, sp, sc)
@testset "Metacommunity" begin
    @test gettypes(meta) == ms
    @test getpartition(meta) == oc_count
    @test isnull(meta.ordinariness)
    @test getordinariness!(meta) ≈ [0.3, 0.6, 1.0]
    @test !isnull(meta.ordinariness)
    @test_throws TypeError Metacommunity(ab3, ms, sc)
    @test_throws ErrorException Metacommunity(-abf, ms, sc)
    @test getsimilarity(gettypes(meta2)) ≈ eye(size(ab3, 1))
    @test Diversity.floattypes(meta) ⊆ mapreduce(Diversity.floattypes, ∩,
                                                 [getabundance(meta),
                                                  getpartition(meta),
                                                  gettypes(meta)])
end

end

