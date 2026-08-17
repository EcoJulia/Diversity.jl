# SPDX-License-Identifier: BSD-2-Clause

module TestMetacommunity
using Test
using LinearAlgebra

using Diversity
using Diversity.API
using Missings

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
    @test floattypes(meta) ⊆ mapreduce(floattypes, ∩,
                    [getabundance(meta),
                        getpartition(meta),
                        gettypes(meta)])
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

end
