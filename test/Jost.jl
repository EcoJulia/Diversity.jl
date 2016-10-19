module TestJost

using Diversity
using Diversity.ᾱ, Diversity.ρ̄

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

# Checking Jost's diversities
using Diversity.Jost

qs = [0, 1, 2, 3, 4, 5, 6, Inf];

numspecies = 100;
numcommunities = 8;
communities = rand(numspecies, numcommunities);
communities /= sum(communities);
probs = reshape(mapslices(sum, communities, 2), (size(communities, 1)));
colweights = rand(numcommunities);
colweights /= sum(colweights);
allthesame = probs * colweights';

@testset "Jost" begin
    @test jostbeta(communities, 1) ≈ (1.0 ./
                                      metadiv(ρ̄(Metacommunity(communities)),
                                               1))
    @test jostbeta(allthesame, qs) ≈ ones(qs)
    
    ## Check Jost's alpha diversity works for all the same subcommunity
    @test jostalpha(allthesame, qs) ≈ metadiv(ᾱ(Metacommunity(allthesame)),
                                               qs)
    
    ## And for all different subcommunities and any subcommunities with the same sizes
    weights = rand(numspecies);
    weights /= sum(weights);
    communitylist = rand(1:numcommunities, numspecies)
    distinct = zeros(Float64, (numspecies, numcommunities))
    for i in 1:numspecies
        distinct[i, communitylist[i]] = weights[i]
    end
    evendistinct = mapslices((x) -> x / (sum(x) * numcommunities), distinct, 1)
    
    @test jostalpha(evendistinct, qs) ≈ metadiv(ᾱ(Metacommunity(evendistinct)), qs)
    
    # Now some even communities, should see that raw and normalised
    # diversities are the same
    smoothed = communities ./ mapslices(sum, communities, 1);
    smoothed /= numcommunities;
    @test jostalpha(smoothed, qs) ≈ metadiv(ᾱ(Metacommunity(smoothed)), qs)
end

end
