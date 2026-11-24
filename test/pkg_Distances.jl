# SPDX-License-Identifier: BSD-2-Clause

module ValidateDistances
using Test
using Statistics

using Diversity
using Diversity.ShortNames
using Diversity.Ecology: pairwisedissimilarity
using Distances

@testset "RenyiDivergence" begin
    @testset "Random RenyiDivergence $i" for i in 1:20
        types = rand(1:(i * 10))
        sc = rand(1:(i * 10))
        pops = rand(types, sc)
        # Make sure not to remove all of the non-zeros from any column
        for j in 1:sc
            vals = pops[:, j] .< median(pops[:, j]) / 2
            for k in 1:types
                if vals[k]
                    pops[k, j] = 0.0
                end
            end
        end
        pops /= sum(pops)
        meta = Metacommunity(pops)
        metapop = getmetaordinariness!(meta)
        weights = getweight(meta)
        pops_norm = pops ./ weights'
        nb = β̄(meta)
        nr = ρ̄(meta)
        b = β(meta)
        r = ρ(meta)
        for q in [rand(7) * 10..., 0, 1, Inf]
            @test exp.(colwise(RenyiDivergence(q),
                               pops_norm, metapop)) ≈
                  subdiv(nb, q)[!, :diversity]
            @test exp.(-colwise(RenyiDivergence(q),
                                pops_norm, metapop)) ≈
                  subdiv(nr, q)[!, :diversity]
            @test exp.(colwise(RenyiDivergence(q),
                               pops_norm, metapop)) ≈
                  subdiv(b, q)[!, :diversity] ./ weights
            @test exp.(-colwise(RenyiDivergence(q),
                                pops_norm, metapop)) ≈
                  subdiv(r, q)[!, :diversity] .* weights
        end
    end
end

# pairwisedissimilarity() against Distances.jl. Every form is a function of
# G = Pᵀ Z P. Writing Z = F Fᵀ, the cosines between nonnegative feature vectors
# of the types (positive semi-definite, in [0, 1], ones on the diagonal), makes
# G the inner products of the subcommunities' profiles Fᵀ P, which Distances.jl
# can measure. Ochiai is their cosine distance, and Rao is half the squared
# Mahalanobis distance in Z between relative abundances. With Z = I and
# presence/absence data, Jaccard and Sørensen are the classic binary indices,
# which are Distances.jl's Jaccard and Bray-Curtis.
@testset "pairwisedissimilarity" begin
    @testset "Random pairwisedissimilarity $i" for i in 1:20
        types = rand(2:(i * 5))
        sc = rand(2:(i * 5))
        ab = rand(types, sc)
        # Make sure not to remove all of the non-zeros from any column
        for j in 1:sc
            ab[ab[:, j] .< median(ab[:, j]) / 2, j] .= 0.0
        end
        ab ./= sum(ab)
        p = ab ./ sum(ab, dims = 1)
        # A small feature of its own keeps every type distinct
        nfeatures = rand(1:(2 * types))
        F = hcat(0.1 .* ((1:types) .== (1:types)'),
                 rand(types, nfeatures) .* (rand(types, nfeatures) .< 0.5))
        F ./= sqrt.(sum(abs2, F, dims = 2))
        Z = F * F'
        Z = (Z + Z') / 2
        for k in 1:types
            Z[k, k] = 1.0
        end

        # Z is positive semi-definite by construction, but often singular, so
        # rounding can fail Distances.jl's check
        sqmahal = pairwise(SqMahalanobis(Z; skipchecks = true), p, dims = 2)

        # Cosines ignore the scale of each subcommunity, so Ochiai is the same
        # with relative abundances
        cosines = pairwise(CosineDist(), F' * ab, dims = 2)
        d = pairwisedissimilarity(ab, Z)
        @test d.ochiai ≈ cosines atol = 1e-12
        @test pairwisedissimilarity(ab, Z; relative = true).ochiai ≈
              cosines atol = 1e-12
        @test d.rao ≈ sqmahal ./ 2 atol = 1e-12
        @test pairwisedissimilarity(ab, Z; squareroot = true).rao ≈
              sqrt.(max.(sqmahal ./ 2, 0)) atol = 1e-6

        # No similarity between types
        @test pairwisedissimilarity(ab).ochiai ≈
              pairwise(CosineDist(), ab, dims = 2) atol = 1e-12
        @test pairwisedissimilarity(ab).rao ≈
              pairwise(SqEuclidean(), p, dims = 2) ./ 2 atol = 1e-12
        pa = Float64.(ab .> 0)
        d01 = pairwisedissimilarity(pa ./ sum(pa))
        @test d01.jaccard ≈ pairwise(Jaccard(), pa, dims = 2) atol = 1e-12
        @test d01.sorensen ≈ pairwise(BrayCurtis(), pa, dims = 2) atol = 1e-12
        @test d01.ochiai ≈ pairwise(CosineDist(), pa, dims = 2) atol = 1e-12
    end
end

end
