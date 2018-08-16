module TestMetacommunity
using Compat.Test
using Compat.LinearAlgebra
using Compat

using Diversity
using Diversity.API
using DataFrames
using Missings

three = [0.3
         0.3
         0.4]
oc_count = Onecommunity()
ab3 = [1 3 0
       2 0 4]
abnorm = ab3 / sum(ab3)
sc = Subcommunities(size(ab3, 2))
@testset "Communities" begin
    oc_2 = Onecommunity("All of it")
    @test countsubcommunities(Onecommunity("Hello")) == countsubcommunities(oc_2)
    @test getsubcommunitynames(oc_2) == ["All of it"]
    @test getsubcommunitynames(oc_count) == ["1"]
    @test getsubcommunitynames(sc) == map(x -> "$x", 1:countsubcommunities(sc))
    @test getsubcommunitynames(Subcommunities(["a", "b"])) == ["a", "b"]
end

sim = [1.0 0.0 0.0
       1.0 1.0 0.0
       1.0 1.0 1.0]
ms = GeneralTypes(sim)
@testset "GeneralTypes" begin
    @test calcsimilarity(ms, 1.0) == sim
    @test_throws DomainError GeneralTypes(-sim)
    @test_throws MethodError GeneralTypes(ab3)
    @test_throws DimensionMismatch GeneralTypes(convert(Matrix{Float64}, ab3))
end

tax = Taxonomy(DataFrame(Species=["This", "That"]), Dict(:Species=>1.0))
@testset "Taxonomy" begin
    @test counttypes(tax) == 2
    @test gettypenames(tax) == ["This", "That"]
    @test_throws ErrorException calcsimilarity(tax, 1.0)
    @test_throws ErrorException _calcordinariness(tax, abnorm, 1.0)
    @test_throws ErrorException Taxonomy(DataFrame(Species=["This", "That"]),
                                         Dict(:Species=>1.0, :Genus=>0.5))
    @test_throws ErrorException Taxonomy(DataFrame(Genus=["This", "That"]),
                                         Dict(:Genus=>1.0))
    @test gettypenames(Taxonomy(DataFrame(Genus=["This", "That"]),
                                Dict(:Genus=>1.0), :Genus)) == ["This", "That"]
end

@testset "Type names" begin
    @test gettypenames(Species(3)) == map(x -> "$x", 1:3)
    @test gettypenames(Species(["My species"])) == ["My species"]
    @test gettypenames(GeneralTypes(Diagonal([1.0]), ["My species"])) == ["My species"]
    @test gettypenames(UniqueTypes(["One", "Two"])) == ["One", "Two"]
end

meta = Metacommunity(three, ms, oc_count)
sp = Species(size(ab3, 1))
abf = ab3 ./ sum(ab3)
meta2 = Metacommunity(abf, sp, sc)
g2 = GeneralTypes(Matrix(1.0I, 2, 2))
@testset "Metacommunity" begin
    @test gettypes(meta) == ms
    @test getpartition(meta) == oc_count
    @test ismissing(meta.ordinariness)
    @test getordinariness!(meta) ≈ [0.3, 0.6, 1.0]
    @test !ismissing(meta.ordinariness)
    @test getabundance(Metacommunity(ab3, g2, sc)) ≈
        getabundance(Metacommunity(abf, g2, sc))
    if VERSION < v"0.7.0-"
        @test_warn "not normalised" Metacommunity(ab3 * 1.0, Matrix(1.0I, 2, 2))
        @test_warn "not normalised" Metacommunity([1.0, 2.0], Matrix(1.0I, 2, 2))
        @test_warn "not normalised" Metacommunity(ab3 * 1.0, meta2)
    end
    @test_nowarn Metacommunity([0.5, 0.5], Matrix(1.0I, 2, 2))
    @test_throws ErrorException Metacommunity(abf, ms, sc)
    @test_throws DimensionMismatch Metacommunity([1, 2, 3]/6, meta2)
    @test_throws DimensionMismatch getabundance(Metacommunity([0.5, 0.5], meta2))
    @test getabundance(Metacommunity(abf, Matrix(1.0I, 2, 2))) ≈
        getabundance(Metacommunity(abf, meta2))
    #@test_throws ErrorException Metacommunity(-abf, g2, sc)
    @test calcsimilarity(gettypes(meta2), _getscale(meta2)) ≈ Matrix(1.0I, size(ab3, 1), size(ab3, 1))
    @test Diversity.floattypes(meta) ⊆ mapreduce(Diversity.floattypes, ∩,
                                                 [getabundance(meta),
                                                  getpartition(meta),
                                                  gettypes(meta)])
end

end
