module TestPhylogenetics
using Base.Test

using PhyloTrees
using Diversity
using Diversity.Phylogenetics

@testset "PhyloTrees" begin
    species = ["Dog", "Cat", "Human"]
    nt = NodeTree(species, nodetype=Vector{Float64})
    n = addnode!(nt)
    addbranch!(nt, n, "Dog", 1.0)
    addbranch!(nt, n, "Cat", 1.0)
    r = addnode!(nt)
    addbranch!(nt, r, "Human", 2.0)
    addbranch!(nt, r, n, 1.0)
    ph = Phylogeny(nt)
    @test Set(species) == Set(getnames(ph))
end

end
