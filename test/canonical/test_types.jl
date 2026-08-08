# SPDX-License-Identifier: BSD-2-Clause
#
# Canonical results for **type construction**: what each `AbstractTypes` subtype turns its input into
# before any measure sees it. `test_measures.jl` pins the arithmetic downstream of a similarity
# matrix; this file pins the matrices themselves, so when a number moves, which file moved says
# whether to look at the measures or at the type.
#
# ⭐ It is also the only numeric gate on the extensions. `ext_*.jl` checks them against hand-computed
# expectations for one tiny input; nothing else records what they produce.

module CanonicalTypes

using Test
using Diversity
using Diversity.ShortNames
using Phylo
using PopGen
using BioSequences

include("canonical.jl")
using .Canonical

const QS = [0, 1, 2, Inf]
qlabel(q) = isinf(q) ? "qInf" : "q$(Int(q))"

@testset "canonical: types" begin
    # ⚠️ The tree is built explicitly rather than drawn with `rand(Nonultrametric(…))`. A blessed
    # value must be a pure function of the code; a random topology would re-bless to noise every run.
    @testset "PhyloBranches" begin
        species = ["Dog", "Human", "Cat"]
        tree = RootedTree(species)
        internal = createnode!(tree)
        createbranch!(tree, internal, species[1], 1.0)
        createbranch!(tree, internal, species[2], 1.0)
        root = createnode!(tree)
        createbranch!(tree, root, internal, 1.0)
        createbranch!(tree, root, species[3], 2.0)

        ph = PhyloBranches(tree)
        meta = Metacommunity([0.4, 0.3, 0.3], ph)

        # ⚠️ `vec` throughout: a metacommunity built from a vector still stores its abundances as a
        # one-column matrix, and the similarity matrix is 5×5. Shapes are asserted below instead.
        blessed("types/phylo/abundance", vec(getabundance(meta)))
        blessed("types/phylo/ordinariness", vec(getordinariness!(meta)))
        blessed("types/phylo/similarity", vec(calcsimilarity(ph, meta.scale)))
        for q in QS
            blessed("types/phylo/gamma_$(qlabel(q))",
                    metadiv(Γ(meta), q)[1, :diversity])
        end

        # ⭐ Raw types are the leaves, processed types the ancestral branches — the distinction the
        # whole `raw::Bool` argument exists for. Shape is asserted here because the blessed vectors
        # above are flat.
        @test gettypenames(ph, true) == species
        @test counttypes(ph, true) == 3
        @test counttypes(ph, false) == 5
        @test length(calcsimilarity(ph, meta.scale)) == 5 * 5
    end

    @testset "GeneticVCF" begin
        pd = PopGen.vcf(Diversity.path("data", "biallelic.vcf"); silent = true)
        gv = GeneticType(pd; distance = :manhattan)
        meta = Metacommunity([0.2 0.1; 0.1 0.3; 0.2 0.1], gv)

        blessed("types/genetic_vcf/similarity", vec(calcsimilarity(gv, 1.0)))
        for q in QS
            blessed("types/genetic_vcf/gamma_$(qlabel(q))",
                    metadiv(Γ(meta), q)[1, :diversity])
        end

        Z = calcsimilarity(gv, 1.0)
        @test size(Z) == (3, 3)
        @test Z == Z'                            # dosage distance is symmetric
        @test all(≈(1.0), Z[i, i] for i in 1:3)  # every sample is identical to itself
    end

    @testset "GeneticFASTA" begin
        seqs = [dna"ACGTACGT", dna"ACGAACGT", dna"TTTTTTTT"]
        gf = GeneticType(seqs; names = ["a", "b", "c"])
        meta = Metacommunity([0.3, 0.3, 0.4], gf)

        blessed("types/genetic_fasta/similarity", vec(calcsimilarity(gf, 1.0)))
        for q in QS
            blessed("types/genetic_fasta/gamma_$(qlabel(q))",
                    metadiv(Γ(meta), q)[1, :diversity])
        end

        Z = calcsimilarity(gf, 1.0)
        @test size(Z) == (3, 3)
        @test Z == Z'
        @test all(≈(1.0), Z[i, i] for i in 1:3)
    end

    # A plain `GeneralTypes` is the control: no derivation at all, so its similarity matrix must come
    # back exactly as handed in. If this moves, the problem is not in any extension.
    @testset "GeneralTypes" begin
        Z = [1.0 0.5 0.0; 0.2 1.0 0.4; 0.0 0.3 1.0]
        gt = GeneralTypes(Z)
        @test calcsimilarity(gt, 1.0) == Z
        blessed("types/general/similarity", vec(calcsimilarity(gt, 1.0)))
    end
end

end
