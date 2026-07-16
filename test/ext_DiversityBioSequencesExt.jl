# SPDX-License-Identifier: BSD-2-Clause

module TestBioSequences
using Test

# Deliberately NOT loading StringDistances here: it ships with Diversity as a
# hard dependency, so the sequence path needs only Diversity + BioSequences.
using Diversity
using Diversity.ShortNames
using BioSequences

@testset "Genetic (sequence)" begin
    seqs = [dna"ACGTACGT", dna"ACGAACGT", dna"TTTTTTTT"]
    seqnames = ["a", "b", "c"]
    gf = GeneticType(seqs; names = seqnames)
    @test gettypenames(gf, true) == seqnames
    @test counttypes(gf, true) == 3

    # Hamming distances: d(a,b)=1, d(a,c)=6, d(b,c)=7; max 7, linear dist2sim.
    dist = [0 1 6; 1 0 7; 6 7 0]
    @test calcsimilarity(gf, 1.0) ≈ 1.0 .- dist ./ 7

    Z = calcsimilarity(gf, 1.0)
    @test Z == Z'                             # symmetric
    @test all(≈(1.0), Z[i, i] for i in 1:3)   # identical to self

    meta = Metacommunity([0.3, 0.3, 0.4], gf)
    @test all(isfinite, metadiv(Γ(meta), 2)[!, :diversity])
end

end
