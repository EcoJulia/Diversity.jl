module TestEcology
using Test

# Looking at relations to historical measures, updated with similarity
# and partitioning
using Diversity
using Diversity.Ecology
using Diversity.ShortNames
using LinearAlgebra

numspecies = 100;
numcommunities = 8;
communities = rand(numspecies, numcommunities);
communities /= sum(communities);
eco = Metacommunity(communities)
weights = rand(numspecies);
weights /= sum(weights);
Z1 = ones(typeof(weights[1]), (length(weights), length(weights)));

@testset "Standard ecological diversities" begin
    @test richness(communities)[!,:diversity] ≈ subdiv(ᾱ(eco), 0)[!,:diversity]
    @test richness(Metacommunity(communities))[!,:diversity] ==
        richness(communities)[!,:diversity]
    @test_throws ErrorException richness(Metacommunity(communities,
        GeneralTypes(rand(numspecies,numspecies))))

    @test shannon(communities)[!,:diversity] ≈ log.(subdiv(ᾱ(eco), 1)[!,:diversity])
    @test shannon(Metacommunity(communities))[!,:diversity] ==
        shannon(communities)[!,:diversity]
    @test_throws ErrorException shannon(Metacommunity(communities,
        GeneralTypes(rand(numspecies,numspecies))))

    @test simpson(communities)[!,:diversity] ≈ subdiv(ᾱ(eco), 2)[!,:diversity] .^ -1
    @test simpson(Metacommunity(communities))[!,:diversity] ==
        simpson(communities)[!,:diversity]
    @test_throws ErrorException simpson(Metacommunity(communities,
        GeneralTypes(rand(numspecies,numspecies))))

    @test jaccard([1 0 0; 0 1 1]'/3)[!,:diversity][1] ≈ 1.0
    @test jaccard(Metacommunity([1 0 1; 0 1 1]'/4))[!,:diversity] ==
        jaccard([1 0 1; 0 1 1]'/4)[!,:diversity]
    @test_throws ErrorException jaccard(Metacommunity([1 0 1; 0 1 1]'/4,
                                                      GeneralTypes(rand(3, 3))))
    @test jaccard([1 0 1; 0 1 1]'/4)[!,:diversity][1] ≈ 2.0 / 3.0
    @test_throws ErrorException jaccard([1 1 0; 0 1 1; 1 1 1]/7)
    @test_throws ErrorException jaccard(Metacommunity([1 1 0; 0 1 1; 1 1 1]/7))

    @test pielou([1, 1])[!,:diversity][1] ≈ 1.0
    @test all(pielou([1 2; 1 2])[!,:diversity] .≈ 1.0)
    @test all(pielou(Metacommunity([1 2; 1 2]))[!,:diversity] .≈ 1.0)

    @test_throws "Can only calculate" gower([1 1 0; 0 1 1; 1 1 1])
    @test gower([1 0; 0 1]).diversity[1] ≈ 1.0
    @test gower([1 1; 1 1]).diversity[1] ≈ 0.0
    @test gower([1 0; 0 0; 0 1; 0 0], countzeros = true).diversity[1] ≈ 0.5
    @test gower([1 0; 0 0; 0 1; 0 0], countzeros = false).diversity[1] ≈ 1.0
end

@testset "Generalised ecological diversities" begin
    @test generalisedrichness(metacommunityDiversity,
                              communities, Z1)[!,:diversity] ≈ [1]
    @test generalisedrichness(subcommunityDiversity,
                              communities)[!,:diversity] ≈
          richness(communities)[!,:diversity]

    @test_throws ErrorException generalisedrichness(individualDiversity,
                                                    communities, Z1)

    @test generalisedshannon(metacommunityDiversity,
                             communities, Z1)[!,:diversity] .+ 1.0 ≈ [1.0]
    @test generalisedshannon(subcommunityDiversity,
                             communities)[!,:diversity] ≈
          shannon(communities)[!,:diversity]
    @test_throws ErrorException generalisedshannon(individualDiversity,
                                                   communities, Z1)

    @test generalisedsimpson(metacommunityDiversity,
                             communities, Z1)[!,:diversity] ≈ [1.0]
    @test generalisedsimpson(subcommunityDiversity,
                             communities)[!,:diversity] ≈
          simpson(communities)[!,:diversity]
    @test_throws ErrorException generalisedsimpson(individualDiversity,
                                                   communities, Z1)

    @test generalisedjaccard([1 0 1; 0 1 1]'/4, UniqueTypes(3))[!,:diversity][1] ≈ 2.0/3.0
    @test generalisedjaccard([1 0 1; 0 1 1]'/4, ones(Float64, 3, 3))[!,:diversity][1] ≈ 0.0
    @test generalisedjaccard([1 1 1; 1 1 1]'/6, Matrix(1.0I, 3, 3))[!,:diversity][1] ≈ 0.0

    @test all(Diversity.Ecology.generalisedpielou(subcommunityDiversity, [1/6 2/6; 1/6 2/6]).diversity .≈ 1.0)

    mat = reshape(rand(9), 3, 3)
    mat ./= sum(mat)
    mats = Float64[1 0 3; 3 2 0; 4 2 3] ./ 18
    matm = Float64[1 0 1; 1 2 0; 2 2 1] ./ 10
    sim = Float64[1 1 0; 1 1 0; 0 0 1]
    @test all(Diversity.Ecology.generalisedpielou(subcommunityDiversity, mat).diversity .≈ pielou(mat).diversity)
    @test all(Diversity.Ecology.generalisedpielou(subcommunityDiversity, mat, UniqueTypes(3)).diversity .≈ pielou(mat).diversity)
    @test_broken all(Diversity.Ecology.generalisedpielou(subcommunityDiversity, mats, sim).diversity .≈ 1.0)
    @test_broken all(Diversity.Ecology.generalisedpielou(metacommunityDiversity, matm, sim).diversity .≈ 1.0)
    @test_throws "Can't calculate Pielou" Diversity.Ecology.generalisedpielou(individualDiversity, mat, sim)
    @test_throws "function cannot run with" pielou(Metacommunity(mat, GeneralTypes(sim)))
end

end
