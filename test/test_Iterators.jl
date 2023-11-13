module TestIterators
using Test
using Phylo
using Diversity

@testset "Iterators" begin
    species = 10
    sc = 5
    abund = rand(species, sc)
    abund ./= sum(abund)
    ru = rand(Ultrametric(species))
    m = Metacommunity(abund, PhyloBranches(ru))
    ti2 = TypeIterator(m)
    ti = TypeIterator(getmetaabundance, m)
    @test length(ti2) == counttypes(PhyloBranches(ru))
    @test length(ti) == counttypes(PhyloBranches(ru))
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
end

end
