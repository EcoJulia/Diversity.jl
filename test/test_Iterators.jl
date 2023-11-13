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
    @test length(TypeIterator(m)) == counttypes(PhyloBranches(ru))
    @test length(TypeIterator(getmetaabundance, m)) == counttypes(PhyloBranches(ru))
    @test length(SubcommunityIterator(m)) == sc
    @test all(reduce(+, SubcommunityIterator(m)) .≈ getmetaabundance(m))
    @test all(reduce(+, TypeIterator(m)) .≈ getweight(m))
    @test all(reduce(+, TypeIterator(getmetaabundance, m))[1] ≈ 1.0)
end

end
