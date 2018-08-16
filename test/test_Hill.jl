module TestHill
using Compat.Test
using Compat.LinearAlgebra

# Checking Hill numbers
using Diversity
using Diversity.Hill

numspecies = 100;
numcommunities = 8;
manyweights = rand(numspecies, numcommunities);
manyweights *= Diagonal(reshape(mapslices(v -> 1. / sum(v), manyweights, dims=1),
                             (numcommunities)));

@testset "Hill numbers" begin
    for i in 1:size(manyweights, 2)
        @test hillnumber(manyweights[:,i], [0])[:diversity] ≈ [numspecies]
    end
end

end
