module TestTypes
using Test
using Diversity
using Diversity.API
using DataFrames
using LinearAlgebra

ab3 = [1 3 0
       2 0 4]
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
    @test getdiversityname(tax) == "Taxonomy"
    @test Float64 ∈ floattypes(tax)
    @test_throws ErrorException calcsimilarity(tax, 1.0)
    # @test_throws ErrorException _calcordinariness(tax, abnorm, 1.0)
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

end
