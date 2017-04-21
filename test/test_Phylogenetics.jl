module TestPhylogenetics
using Base.Test

using PhyloTrees
using Diversity
using Diversity.Phylogenetics

@testset "PhyloTrees" begin
    species = ["Dog", "Human", "Cat"]
    abund = [0.4, 0.3, 0.3]
    nt = NamedTree(species)
    n = addnode!(nt)
    addbranch!(nt, n, "Dog", 1.0)
    addbranch!(nt, n, "Cat", 1.0)
    r = addnode!(nt)
    addbranch!(nt, r, "Human", 2.0)
    addbranch!(nt, r, n, 1.0)
    ph = Phylogeny(nt)
    leafnames = getphylonames(ph)
    order = mapreduce(name -> find(species .== name), append!, leafnames)
    @test species[order] == getphylonames(ph)
    metaphylo = Metacommunity(abund[order], ph)
    @test getabundance(metaphylo) ≈ [0.15, 0.2, 0.3, 0.15, 0.2]
    @test getordinariness!(metaphylo) ≈ [0.7, 0.7, 0.3, 0.3, 0.4]
    @test meta_gamma(metaphylo, 0)[:diversity] == [2.5]
    @test sub_gamma(metaphylo, 0)[:diversity] == [2.5]
end

end
