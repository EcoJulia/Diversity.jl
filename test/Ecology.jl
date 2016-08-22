module TestEcology

using Diversity
using Diversity.ᾱ, Diversity.γ
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

# Looking at relations to historical measures, updated with similarity
# and partitioning
using Diversity.Ecology

numspecies = 100;
numcommunities = 8;
communities = rand(numspecies, numcommunities);
communities /= sum(communities);
eco = Supercommunity(communities)
weights = rand(numspecies);
weights /= sum(weights);
Z1 = ones(typeof(weights[1]), (length(weights), length(weights)));

@testset "Standard ecological diversities" begin
    @test richness(communities) ≈ subdiv(ᾱ(eco), 0)
    
    @test shannon(communities) ≈ log(subdiv(ᾱ(eco), 1))
    
    @test simpson(communities) ≈ subdiv(ᾱ(eco), 2) .^ -1
    
    @test jaccard([1 0 0; 0 1 1]') ≈ 0
    @test jaccard([1 0 1; 0 1 1]') ≈ 1.0 / 3.0
    @test_throws ErrorException jaccard([1 1 0; 0 1 1; 1 1 1])
end

@testset "Generalised ecological diversities" begin
    @test generalisedrichness(supercommunityDiversity, γ,
                              communities, Z1) ≈ 1
    
    @test generalisedshannon(supercommunityDiversity, γ,
                             communities, Z1) ≈ 0
    
    @test generalisedsimpson(supercommunityDiversity, γ,
                             communities, Z1) ≈ 1
    
    @test generalisedjaccard([1 0 1; 0 1 1]', [0, Inf]) ≈ [1.0/3.0, 1.0]
    @test generalisedjaccard([1 1 1; 1 1 1]', [0, 1]) ≈ [1, 1]
end

end
