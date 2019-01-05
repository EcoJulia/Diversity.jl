module TestJost
using Compat.Test
using Compat

# Checking Jost's diversities
using Diversity
using Diversity.ShortNames
using Diversity.Jost

qs = [0, 1, 2, 3, 4, 5, 6, Inf];

numspecies = 100;
numcommunities = 8;
communities = rand(numspecies, numcommunities);
communities /= sum(communities);
probs = reshape(mapslices(sum, communities, dims=2), (size(communities, 1)));
colweights = rand(numcommunities);
colweights /= sum(colweights);
allthesame = probs * colweights';

@testset "Jost" begin
    @test jostbeta(communities, 1)[:diversity] ≈
        (1.0 ./ metadiv(ρ̄(Metacommunity(communities)), 1)[:diversity])
    @test jostbeta(allthesame, qs)[:diversity] ≈ fill!(similar(qs), 1)
    @test jostbeta(Metacommunity(allthesame), 1)[:diversity] ==
        jostbeta(allthesame, 1)[:diversity]
    @test_throws ErrorException jostbeta(Metacommunity(allthesame,
        GeneralTypes(rand(numspecies, numspecies))), 1)
    ## Check Jost's alpha diversity works for all the same subcommunity
    @test jostalpha(allthesame, qs)[:diversity] ≈
          metadiv(ᾱ(Metacommunity(allthesame)), qs)[:diversity]
    @test jostalpha(Metacommunity(allthesame), 1)[:diversity] ==
        jostalpha(allthesame, 1)[:diversity]
    @test_throws ErrorException jostalpha(Metacommunity(allthesame,
        GeneralTypes(rand(numspecies, numspecies))), 1)

    ## And for all different subcommunities and any subcommunities with the same sizes
    weights = rand(numspecies);
    weights /= sum(weights);
    communitylist = rand(1:numcommunities, numspecies)
    distinct = zeros(Float64, (numspecies, numcommunities))
    for i in 1:numspecies
        distinct[i, communitylist[i]] = weights[i]
    end
    evendistinct = mapslices((x) -> x / (sum(x) * numcommunities), distinct, dims=1)

    @test jostalpha(evendistinct, qs)[:diversity] ≈
          metadiv(ᾱ(Metacommunity(evendistinct)), qs)[:diversity]

    # Now some even communities, should see that raw and normalised
    # diversities are the same
    smoothed = communities ./ mapslices(sum, communities, dims=1);
    smoothed /= numcommunities;
    @test jostalpha(smoothed, qs)[:diversity] ≈
          metadiv(ᾱ(Metacommunity(smoothed)), qs)[:diversity]
end

end
