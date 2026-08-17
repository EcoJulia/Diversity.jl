# SPDX-License-Identifier: BSD-2-Clause

module TestIterators
using Test
using Diversity

# Deliberately built from `GeneralTypes` rather than a phylogeny, so this file needs nothing
# beyond the package's own dependencies and runs as a bare script. The iterators do not care what
# kind of types they are walking. The phylogenetic case — where the processed types outnumber the
# raw ones — is covered in `ext_DiversityPhyloExt.jl`, where `Phylo` is loaded anyway.
@testset "Iterators" begin
    species = 10
    sc = 5
    abund = rand(species, sc)
    abund ./= sum(abund)
    Z = fill(0.5, species, species)
    for i in 1:species
        Z[i, i] = 1.0
    end
    m = Metacommunity(abund, GeneralTypes(Z))

    ti2 = TypeIterator(m)
    ti = TypeIterator(getmetaabundance, m)
    @test length(ti2) == counttypes(m)
    @test length(ti) == counttypes(m)
    @test_throws "Can't iterate" TypeIterator(sum ∘ getweight, m)

    si = SubcommunityIterator(m)
    @test length(si) == sc
    @test all(reduce(+, si) .≈ getmetaabundance(m))
    @test all(reduce(+, ti2) .≈ getweight(m))
    @test all(reduce(+, ti)[1] ≈ 1.0)
    @test sum(SubcommunityIterator(getweight, m))[1] ≈ 1.0
    @test_throws "Can't iterate" SubcommunityIterator(sum ∘ getabundance, m)

    @test Base.IteratorSize(typeof(ti)) == Base.HasLength()
    @test Base.IteratorSize(typeof(ti2)) == Base.HasLength()
    @test Base.IteratorEltype(typeof(ti)) == Base.HasEltype()
    @test eltype(ti) ≡ Float64

    @test Base.IteratorSize(typeof(si)) == Base.HasLength()
    @test Base.IteratorEltype(typeof(si)) == Base.HasEltype()
    @test eltype(si) ≡ Float64

    # Iterating types and iterating subcommunities partition the same abundances two ways, so both
    # must total the same thing.
    @test sum(sum, ti2) ≈ sum(sum, si) ≈ 1.0
end

end
