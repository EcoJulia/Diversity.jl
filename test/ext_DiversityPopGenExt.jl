# SPDX-License-Identifier: BSD-2-Clause

module TestPopGen
using Test

using Diversity
using Diversity.ShortNames
using PopGen
using DataFrames

@testset "Genetic (VCF)" begin
    # test/data/biallelic.vcf matches the rdiversity gen2dist test data:
    #   s1 = 0|0, 0|1, 1|1   s2 = 1|1, 0|1, 0|0   s3 = 0|0, 0|0, 0|0
    # gen2dist(biallelic = TRUE) gives Manhattan distances [[0,4,3],[4,0,3],[3,3,0]].
    file = joinpath(@__DIR__, "data", "biallelic.vcf")
    pd = PopGen.vcf(file; silent = true)

    gv = GeneticType(pd; distance = :manhattan)
    @test gettypenames(gv, true) == ["s1", "s2", "s3"]
    @test counttypes(gv, true) == 3

    # Manhattan distances [[0,4,3],[4,0,3],[3,3,0]] normalised by max 4, then
    # linear dist2sim gives Z = 1 - d/4.
    expected = 1.0 .- [0 4 3; 4 0 3; 3 3 0] ./ 4
    @test calcsimilarity(gv, 1.0) ≈ expected

    # Exponential transform matches rdiversity dist2sim(transform = "exponential").
    ge = GeneticType(pd; distance = :manhattan, transform = :exponential)
    @test calcsimilarity(ge, 1.0) ≈ exp.(-([0 4 3; 4 0 3; 3 3 0] ./ 4))

    # vcf_dataframe round-trips the PopData back to VCF-body genotype strings
    # (the structure rdiversity's gen2dist consumes).
    vcf = vcf_dataframe(pd)
    @test names(vcf)[9:end] == ["FORMAT", "s1", "s2", "s3"]
    @test all(==("GT"), vcf.FORMAT)
    @test vcf.s1 == ["0|0", "0|1", "1|1"]
    @test vcf.s2 == ["1|1", "0|1", "0|0"]
    @test vcf.s3 == ["0|0", "0|0", "0|0"]

    # End-to-end through a metacommunity.
    pops = [0.2 0.1; 0.1 0.3; 0.2 0.1]
    meta = Metacommunity(pops, gv)
    @test gettypenames(meta, true) == ["s1", "s2", "s3"]
    @test all(isfinite, metadiv(Γ(meta), 1)[!, :diversity])
    @test all(isfinite, subdiv(ᾱ(meta), 0)[!, :diversity])
end

end
