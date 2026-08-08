# SPDX-License-Identifier: BSD-2-Clause

module TestAPI
using Test

using Diversity
using Diversity.API
using Diversity.ShortNames
using LinearAlgebra

# ⭐ A minimal third-party `AbstractTypes`, implementing *only* the two methods the API documents as
# required. Everything else must come from the defaults — which is the whole claim `Diversity.API`
# makes, and which nothing else in this suite checks.
#
# This matters beyond tidiness: EcoSISTEM attaches its `SpeciesList`, `AbstractHabitat` and
# `Ecosystem` to this same API and is never loaded by this suite, so a change that quietly adds a
# required method would go unnoticed here and break there. This type is the canary.
struct UniformSimilarity <: Diversity.API.AbstractTypes
    names::Vector{String}
    similarity::Float64
end

Diversity.API._gettypenames(us::UniformSimilarity, ::Bool) = us.names

function Diversity.API._calcsimilarity(us::UniformSimilarity, ::Real)
    n = length(us.names)
    Z = fill(us.similarity, n, n)
    for i in 1:n
        Z[i, i] = 1.0
    end
    return Z
end

@testset "floattypes" begin
    @test Float32 ∈ floattypes(Float32[1.0])
    @test Float64 ∈ floattypes(Float64[1.0])

    # An AbstractTypes or AbstractPartition with no opinion accepts every float type; a
    # metacommunity is pinned to the one it was built with.
    @test Float32 ∈ floattypes(UniqueTypes(3))
    @test Float64 ∈ floattypes(UniqueTypes(3))
    @test Float32 ∈ floattypes(Subcommunities(2))
    @test floattypes(Metacommunity(Float64[0.5, 0.5])) == Set([Float64])
    @test floattypes(GeneralTypes(Matrix(1.0I, 2, 2))) == Set([Float64])
end

@testset "typematch" begin
    @test typematch(Float64[1.0], UniqueTypes(1), Onecommunity())
    @test typematch(Float64[1.0], Float64[2.0])
    # Two concrete float types with nothing in common must not match.
    @test !typematch(Float32[1.0], Float64[1.0])
end

@testset "mcmatch" begin
    types = UniqueTypes(3)
    part = Subcommunities(2)
    good = [0.2 0.1; 0.1 0.3; 0.2 0.1]
    @test mcmatch(good, types, part)

    # 🔴 Each failure mode separately — this is the gate every `Metacommunity` constructor runs, and
    # a silent hole in it would let a mis-shaped metacommunity through to the measures.
    @test !mcmatch(good, UniqueTypes(2), part)          # wrong number of types
    @test !mcmatch(good, types, Subcommunities(3))      # wrong number of subcommunities
    @test !mcmatch(good ./ 2, types, part)              # abundances do not sum to 1
    # Float mismatch: `GeneralTypes` is pinned to the element type of its Z, so Float32 abundances
    # against a Float64 similarity matrix have no float type in common.
    @test !mcmatch(Float32.(good), GeneralTypes(Matrix(1.0I, 3, 3)), part)
end

@testset "Minimal AbstractTypes implementation" begin
    types = UniformSimilarity(["a", "b", "c"], 0.5)

    # Not implemented here — supplied by the defaults, which is the point.
    @test counttypes(types) == 3
    @test gettypenames(types) == ["a", "b", "c"]
    # A type that does not name itself is "unknown", not "species": `AbstractTypes` overrides the
    # looser `AbstractThings` default precisely so an un-named type says so in the output DataFrame.
    @test getdiversityname(types) == "unknown"
    @test calcsimilarity(types, 1.0) == [1.0 0.5 0.5; 0.5 1.0 0.5; 0.5 0.5 1.0]

    # It works everywhere a built-in type does.
    meta = Metacommunity([0.5, 0.3, 0.2], types)
    @test counttypes(meta) == 3
    @test countsubcommunities(meta) == 1
    @test getordinariness!(meta) ≈
          [1.0 0.5 0.5; 0.5 1.0 0.5; 0.5 0.5 1.0] *
          reshape([0.5, 0.3, 0.2], 3, 1)
    @test all(isfinite, metadiv(Γ(meta), [0, 1, 2, Inf])[!, :diversity])

    # ⭐ The similarity parameter brackets the answer: at 0 the three types are wholly distinct and
    # the metacommunity holds three types' worth of diversity, at 1 they are interchangeable and it
    # holds one. Anything in between must fall between.
    distinct = Metacommunity([0.5, 0.3, 0.2],
                             UniformSimilarity(["a", "b", "c"], 0.0))
    identical = Metacommunity([0.5, 0.3, 0.2],
                              UniformSimilarity(["a", "b", "c"], 1.0))
    @test metadiv(Γ(distinct), 0)[1, :diversity] ≈ 3
    @test metadiv(Γ(identical), 0)[1, :diversity] ≈ 1
    @test 1 < metadiv(Γ(meta), 0)[1, :diversity] < 3
end

end
