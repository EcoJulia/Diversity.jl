module TestPhylogenetics
using Test

using Phylo
using Diversity

@testset "Phylo" begin
    species = ["Dog", "Human", "Cat"]
    abund = [0.4, 0.3, 0.3]
    nt = RootedTree(species)
    n = createnode!(nt)
    createbranch!(nt, n, species[1], 1.0)
    createbranch!(nt, n, species[2], 1.0)
    r = createnode!(nt)
    createbranch!(nt, r, n, 1.0)
    createbranch!(nt, r, species[3], 2.0)
    ph = PhyloBranches(nt)
    leafnames = gettypenames(ph, true)
    @test species == gettypenames(ph, true)
    metaphylo = Metacommunity(abund, ph)
    @test gettypenames(metaphylo, true) == species
    @test getabundance(metaphylo, true) ≈ abund
    @test getabundance(metaphylo) ≈ [0.2, 0.2, 0.15, 0.15, 0.3]
    @test getordinariness!(metaphylo) ≈ [0.4, 0.7, 0.3, 0.7, 0.3]
    @test meta_gamma(metaphylo, 0)[!,:diversity] == [2.5]
    @test sub_gamma(metaphylo, 0)[!,:diversity] == [2.5]
end

end
