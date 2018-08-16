module TestPhylogenetics
using Compat.Test

using Phylo
using Diversity

@testset "Phylo" begin
    species = ["Dog", "Human", "Cat"]
    abund = [0.4, 0.3, 0.3]
    nt = NamedTree(species)
    n = addnode!(nt)
    addbranch!(nt, n, "Dog", 1.0)
    addbranch!(nt, n, "Cat", 1.0)
    r = addnode!(nt)
    addbranch!(nt, r, "Human", 2.0)
    addbranch!(nt, r, n, 1.0)
    ph = PhyloTypes(nt)
    leafnames = gettypenames(ph, true)
    @test species == gettypenames(ph, true)
    metaphylo = Metacommunity(abund, ph)
    @test gettypenames(metaphylo, true) == species
    @test getabundance(metaphylo, true) ≈ abund
    @test getabundance(metaphylo) ≈ [0.15, 0.2, 0.3, 0.15, 0.2]
    @test getordinariness!(metaphylo) ≈ [0.7, 0.7, 0.3, 0.3, 0.4]
    @test meta_gamma(metaphylo, 0)[:diversity] == [2.5]
    @test sub_gamma(metaphylo, 0)[:diversity] == [2.5]
end

end
