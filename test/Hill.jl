module TestHill

using Diversity
using Base.Test

# Checking Hill numbers
using Diversity.Hill

numspecies = 100;
numcommunities = 8;
manyweights = rand(numspecies, numcommunities);
manyweights *= diagm(reshape(mapslices(v -> 1. / sum(v), manyweights, 1),
                             (numcommunities)));

for i in 1:size(manyweights, 2)
    @test_approx_eq hillnumber(manyweights[:,i], [0]) numspecies * ones((1, size(manyweights[:,i], 2)))
end

end
