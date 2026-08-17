# SPDX-License-Identifier: BSD-2-Clause

module TestBioSequences
using Test

# Deliberately NOT loading StringDistances here: it ships with Diversity as a
# hard dependency, so the sequence path needs only Diversity + BioSequences.
using Diversity
using Diversity.ShortNames
using BioSequences
using FASTX

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

@testset "Genetic (sequence) from a FASTA file" begin
    # The documented route into this extension is "a vector of aligned `BioSequence`s", but a user
    # gets those from a file. Reading `test/data/sequences.fasta` — the same three sequences as
    # above — checks that the whole path works and that the names come from the file rather than
    # being invented.
    file = Diversity.path("data", "sequences.fasta")
    records = open(FASTAReader, file) do reader
        return collect(reader)
    end
    seqs = [sequence(LongDNA{4}, record) for record in records]
    gf = GeneticType(seqs; names = identifier.(records))

    @test gettypenames(gf, true) == ["a", "b", "c"]
    @test calcsimilarity(gf, 1.0) ≈ 1.0 .- [0 1 6; 1 0 7; 6 7 0] ./ 7

    # Identical to building the same sequences from literals, which is the claim worth making.
    literal = GeneticType([dna"ACGTACGT", dna"ACGAACGT", dna"TTTTTTTT"];
                          names = ["a", "b", "c"])
    @test calcsimilarity(gf, 1.0) == calcsimilarity(literal, 1.0)
end

end
