module TestEcology

using Diversity
using Base.Test

# Looking at relations to historical measures, updated with similarity
# and partitioning
using Diversity.Ecology

numspecies = 100;
numcommunities = 8;
communities = rand(numspecies, numcommunities);
communities /= sum(communities);
weights = rand(numspecies);
weights /= sum(weights);
Z1 = ones(typeof(weights[1]), (length(weights), length(weights)));

@test_approx_eq richness(communities) Dᾱ(communities, 0)
@test_approx_eq generalisedrichness(DḠ, communities, Z1) 1

@test_approx_eq shannon(communities) log(Dᾱ(communities, 1))
@test_approx_eq generalisedshannon(DḠ, communities, Z1) 0

@test_approx_eq simpson(communities) Dᾱ(communities, 2) .^ -1
@test_approx_eq generalisedsimpson(DḠ, communities, Z1) 1

@test_approx_eq jaccard([1 0 0; 0 1 1.]') 0
@test_approx_eq jaccard([1 0 1; 0 1 1.]') 1 / 3
@test_throws ErrorException jaccard([1 1 0; 0 1 1; 1 1 1.])

@test_approx_eq generalisedjaccard([1 0 1; 0 1 1.]', [0, Inf]) [1/3, 1]
@test_approx_eq generalisedjaccard([1 1 1; 1 1 1.]', [0, 1]) [1, 1]

end
