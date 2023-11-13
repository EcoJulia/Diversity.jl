module TestHill
using Test
using LinearAlgebra

# Checking Hill numbers
using Diversity
using Diversity.Hill

numspecies = 100;
numcommunities = 8;
manyweights = rand(numspecies, numcommunities);
manyweights *= Diagonal(reshape(mapslices(v -> 1. / sum(v), manyweights, dims=1),
                             (numcommunities)));
meta = Metacommunity(manyweights, 1.0I(100))

@testset "Hill numbers" begin
    for i in axes(manyweights, 2)
        @test hillnumber(manyweights[:,i], [0])[!,:diversity] ≈ [numspecies]
    end
    @test_throws "types as contains similarity" hillnumber(meta, 0)
end

end
