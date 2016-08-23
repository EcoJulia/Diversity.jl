module TestHill

using Diversity
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

# Checking Hill numbers
using Diversity.Hill

numspecies = 100;
numcommunities = 8;
manyweights = rand(numspecies, numcommunities);
manyweights *= diagm(reshape(mapslices(v -> 1. / sum(v), manyweights, 1),
                             (numcommunities)));

@testset "Hill numbers" begin
    for i in 1:size(manyweights, 2)
        @test hillnumber(manyweights[:,i], [0]) ≈ numspecies * ones((1, size(manyweights[:,i], 2)))
    end
end

end
