module TestInterface
using Test

# Checking EcoBase interface
using Diversity
using EcoBase

numspecies = 10;
numcommunities = 8;
manyweights = rand(numspecies, numcommunities);
manyweights /= sum(manyweights);

@testset "Text output" begin
    species = map(n -> "Species $n", 1:numspecies)
    communities = map(n -> "SC $n", 1:numcommunities)
    ut = UniqueTypes(species)
    sc = Subcommunities(communities)
    mc = Metacommunity(manyweights, ut, sc)

    io = IOBuffer()
    show(io, mc)
    @test occursin("measuring", String(take!(io)))
end

end
