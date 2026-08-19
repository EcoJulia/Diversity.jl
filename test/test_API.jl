# SPDX-License-Identifier: BSD-2-Clause

module TestAPI
using Test

using Diversity
using Diversity.API
using Diversity.ShortNames
using LinearAlgebra

# A minimal third-party `AbstractTypes`, implementing *only* the two methods the API documents as
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

    # An AbstractTypes or AbstractPartition with no opinion accepts every float type. It says so
    # with the abstract type standing for all of them, rather than by listing the concrete ones --
    # which used to mean enumerating `subtypes(AbstractFloat)`, and so missing any float that was
    # not a direct subtype. Assert the meaning rather than the representation.
    @test all(F -> any(u -> F <: u, floattypes(UniqueTypes(3))),
              (Float16, Float32, Float64, BigFloat))
    @test all(F -> any(u -> F <: u, floattypes(Subcommunities(2))),
              (Float16, Float32, Float64, BigFloat))
    @test typematch(Float32[1.0], UniqueTypes(1), Subcommunities(1))
    @test typematch(BigFloat[1.0], UniqueTypes(1), Subcommunities(1))
    @test floattypes(Metacommunity(Float64[0.5, 0.5])) == Set{Type}([Float64])
    # GeneralTypes is pinned to the float type of its own similarity matrix, so unlike a bare
    # UniqueTypes it does have an opinion about which float it works with.
    @test floattypes(GeneralTypes(Matrix(1.0I, 2, 2))) == Set{Type}([Float64])
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

    # Note: Each failure mode separately — this is the gate every `Metacommunity` constructor runs, and
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

    # The similarity parameter brackets the answer: at 0 the three types are wholly distinct and
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

# A float type that is not a *direct* subtype of AbstractFloat, which is what the old
# implementation enumerated. Declared at module scope because a type cannot be defined inside a
# testset.
abstract type IndirectFloats <: AbstractFloat end
struct Indirect <: IndirectFloats
    x::Float64
end

@testset "Float type compatibility" begin
    # A type or a partition works with any float, and says so with the abstract type rather than
    # by listing the concrete ones.
    @test floattypes(UniqueTypes(2)) == Set{Type}([AbstractFloat])
    @test floattypes(Subcommunities(2)) == Set{Type}([AbstractFloat])
    # An array or a metacommunity is committed to the one it holds.
    @test floattypes(rand(Float64, 2, 2)) == Set{Type}([Float64])
    @test floattypes(rand(Float32, 2, 2)) == Set{Type}([Float32])
    @test floattypes(Metacommunity(rand(2, 2) ./ 2)) == Set{Type}([Float64])

    # An abstract entry stands for any of its subtypes, so "any float" meets Float64 at Float64
    # rather than at nothing.
    @test typematch(rand(2, 2), UniqueTypes(2), Subcommunities(2))
    @test typematch(rand(Float32, 2, 2), UniqueTypes(2), Subcommunities(2))
    # Two different concrete floats are still incompatible, which is the point of the check.
    @test !typematch(rand(Float32, 2, 2),
                     Metacommunity(rand(Float64, 2, 2) ./ 2))
    @test typematch(rand(Float64, 2, 2),
                    Metacommunity(rand(Float64, 2, 2) ./ 2))

    # The case this used to get wrong: enumerating subtypes(AbstractFloat) lists only the direct
    # ones, so a float defined a level deeper was rejected outright and its metacommunity refused.
    @test Indirect <: AbstractFloat
    @test floattypes(Indirect[]) == Set{Type}([Indirect])
    @test typematch(Indirect[], UniqueTypes(2), Subcommunities(2))
end

end
