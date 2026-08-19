# SPDX-License-Identifier: BSD-2-Clause

module TestMetacommunity
using Test
using LinearAlgebra

using Diversity
using Diversity.API
using EcoBase
using Missings

# An assemblage from outside our hierarchy whose *types* nonetheless carry similarity — the shape
# EcoSISTEM's Ecosystem has. Its places are an EcoBase `AbstractPlaces`, which cannot also be an
# `AbstractPartition`, so this is what reaches the second `_aspartition` method.
struct ForeignPlaces <: EcoBase.AbstractPlaces{Nothing} end
EcoBase.placenames(::ForeignPlaces) = ["west", "east"]

struct ForeignAssemblage{T <: Diversity.AbstractTypes} <:
       EcoBase.AbstractAssemblage{Float64, T, ForeignPlaces}
    types::T
    abundances::Matrix{Float64}
end
EcoBase.things(fa::ForeignAssemblage) = fa.types
EcoBase.places(::ForeignAssemblage) = ForeignPlaces()
EcoBase.occurrences(fa::ForeignAssemblage) = fa.abundances

three = [0.3
         0.3
         0.4]
three_1 = [3
           3
           4]
oc_count = Onecommunity()
sim = [1.0 0 0; 1.0 1.0 0.0; 1.0 1.0 1.0]
ms = GeneralTypes(sim)
meta = Metacommunity(three, ms, oc_count)
ab3 = [1 3 0
       2 0 4]
sp = Species(size(ab3, 1))
abf = ab3 ./ sum(ab3)
sc = Subcommunities(size(ab3, 2))
meta2 = Metacommunity(abf, sp, sc)
g2 = GeneralTypes(Matrix(1.0I, 2, 2))
@testset "Metacommunity" begin
    @test meta_gamma(Metacommunity(three_1, meta), 1).diversity[1] ≈
          meta_gamma(Metacommunity(three, meta), 1).diversity[1]
    @test_warn "Abundances not normalised" meta_gamma(Metacommunity(three_1,
                                                                    meta), 0)
    @test_warn "Abundances not normalised" Metacommunity(ab3, meta2)

    @test gettypes(meta) == ms
    @test getpartition(meta) == oc_count
    @test ismissing(meta.ordinariness)
    @test getordinariness!(meta) ≈ [0.3, 0.6, 1.0]
    @test !ismissing(meta.ordinariness)
    @test getabundance(Metacommunity(ab3, g2, sc)) ≈
          getabundance(Metacommunity(abf, g2, sc))
    @test_nowarn Metacommunity([0.5, 0.5], Matrix(1.0I, 2, 2))
    @test_throws ErrorException Metacommunity(abf, ms, sc)
    @test_throws DimensionMismatch Metacommunity([1, 2, 3] / 6, meta2)
    @test_throws DimensionMismatch getabundance(Metacommunity([0.5, 0.5],
                                                              meta2))
    @test getabundance(Metacommunity(abf, Matrix(1.0I, 2, 2))) ≈
          getabundance(Metacommunity(abf, meta2))
    #@test_throws ErrorException Metacommunity(-abf, g2, sc)
    @test calcsimilarity(gettypes(meta2), _getscale(meta2)) ≈
          Matrix(1.0I, size(ab3, 1), size(ab3, 1))
    # The metacommunity's float type has to be compatible with the parts it was built from. Asked
    # through `typematch`, which is the contract: a literal set intersection will not do it, since
    # a types object or a partition answers "any float" with the abstract type rather than by
    # listing the concrete ones.
    @test typematch(getabundance(meta), getpartition(meta), gettypes(meta))
    @test floattypes(meta) == floattypes(getabundance(meta))
end

@testset "Counts with a similarity matrix" begin
    Z = Matrix(1.0I, 2, 2)
    @test getabundance(Metacommunity(ab3, Z)) ≈
          getabundance(Metacommunity(abf, Z))
    Z3 = Matrix(1.0I, 3, 3)      # `three_1` has three types, so it needs a 3x3 similarity matrix
    @test getabundance(Metacommunity(three_1, Z3)) ≈
          getabundance(Metacommunity(three_1 ./ sum(three_1), Z3))

    # Integer counts normalise silently; floats that miss still warn.
    @test_nowarn Metacommunity(ab3, Z)
    @test_warn "Abundances not normalised" Metacommunity(abf .* 2, Z)

    # Genuinely mismatched float types are still refused — but by `mcmatch`, which says so,
    # rather than by there being no applicable method at all.
    @test_throws ErrorException Metacommunity(Float32.(abf), Z)
end

@testset "Translating an assemblage that has similarity" begin
    # `Metacommunity(::AbstractAssemblage)` must reach a metacommunity that computes
    # similarity-sensitive diversity whatever the source types were, by materialising the
    # similarity into a GeneralTypes.
    Z = [1.0 0.5 0.0; 0.5 1.0 0.5; 0.0 0.5 1.0]
    types = GeneralTypes(Z, ["ash", "oak", "elm"])
    part = Subcommunities(["north", "south"])
    source = Metacommunity([0.1 0.2; 0.2 0.1; 0.2 0.2], types, part)
    conv = Metacommunity(source)

    @test gettypes(conv) isa GeneralTypes
    @test calcsimilarity(gettypes(conv), 1) ≈ Z
    # The names of both components survive the translation.
    @test gettypenames(conv) == gettypenames(source)
    @test getsubcommunitynames(conv) == getsubcommunitynames(source)

    # And it is the *same* metacommunity as far as every measure is concerned.
    for q in [0, 1, 2, Inf]
        @test norm_sub_alpha(conv, q).diversity ≈
              norm_sub_alpha(source, q).diversity
        @test norm_sub_beta(conv, q).diversity ≈
              norm_sub_beta(source, q).diversity
        @test norm_sub_rho(conv, q).diversity ≈
              norm_sub_rho(source, q).diversity
        @test meta_gamma(conv, q).diversity ≈ meta_gamma(source, q).diversity
    end

    # Without similarity the other branch still gives UniqueTypes, as it always has.
    plain = Metacommunity([0.1 0.2; 0.2 0.1; 0.2 0.2])
    @test gettypes(Metacommunity(plain)) isa UniqueTypes
end

@testset "Translating a foreign assemblage that has similarity" begin
    # The other `_aspartition` method: the source is an EcoBase assemblage from outside this
    # package, so its places cannot be an AbstractPartition and a Subcommunities of the right size
    # is built instead. Everything else must still come across.
    Z = [1.0 0.5 0.0; 0.5 1.0 0.5; 0.0 0.5 1.0]
    types = GeneralTypes(Z, ["ash", "oak", "elm"])
    foreign = ForeignAssemblage(types, [0.1 0.2; 0.2 0.1; 0.2 0.2])
    conv = Metacommunity(foreign)

    @test getpartition(foreign) isa ForeignPlaces
    @test !(getpartition(foreign) isa Diversity.AbstractPartition)
    @test getpartition(conv) isa Subcommunities
    @test countsubcommunities(conv) == 2

    # The similarity and both sets of names survive: the partition object itself cannot be reused,
    # but its `placenames` are carried into the Subcommunities built to replace it.
    @test gettypes(conv) isa GeneralTypes
    @test calcsimilarity(gettypes(conv), 1) ≈ Z
    @test gettypenames(conv) == ["ash", "oak", "elm"]
    @test getsubcommunitynames(conv) == ["west", "east"]
    @test getsubcommunitynames(conv) == placenames(getpartition(foreign))

    # And it measures the same as the assemblage it came from, which is the point of translating.
    for q in [0, 1, 2, Inf]
        @test norm_sub_alpha(conv, q).diversity ≈
              norm_sub_alpha(foreign, q).diversity
        @test norm_sub_rho(conv, q).diversity ≈
              norm_sub_rho(foreign, q).diversity
        @test meta_gamma(conv, q).diversity ≈ meta_gamma(foreign, q).diversity
    end
end

end
