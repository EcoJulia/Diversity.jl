# SPDX-License-Identifier: BSD-2-Clause

module TestEcology
using Test

# Looking at relations to historical measures, updated with similarity
# and partitioning
using Diversity
using Diversity.Ecology
using Diversity.Ecology: generalisedpielou
using Diversity.ShortNames
using LinearAlgebra

numspecies = 100
numcommunities = 8
communities = rand(numspecies, numcommunities)
communities /= sum(communities)
eco = Metacommunity(communities)
weights = rand(numspecies)
weights /= sum(weights)
Z1 = ones(typeof(weights[1]), (length(weights), length(weights)))

@testset "Standard ecological diversities" begin
    @test richness(communities)[!, :diversity] ≈
          subdiv(ᾱ(eco), 0)[!, :diversity]
    @test richness(Metacommunity(communities))[!, :diversity] ==
          richness(communities)[!, :diversity]
    @test_throws ErrorException richness(Metacommunity(communities,
                                                       GeneralTypes(rand(numspecies,
                                                                         numspecies))))

    @test shannon(communities)[!, :diversity] ≈
          log.(subdiv(ᾱ(eco), 1)[!, :diversity])
    @test shannon(Metacommunity(communities))[!, :diversity] ==
          shannon(communities)[!, :diversity]
    @test_throws ErrorException shannon(Metacommunity(communities,
                                                      GeneralTypes(rand(numspecies,
                                                                        numspecies))))

    @test simpson(communities)[!, :diversity] ≈
          subdiv(ᾱ(eco), 2)[!, :diversity] .^ -1
    @test simpson(Metacommunity(communities))[!, :diversity] ==
          simpson(communities)[!, :diversity]
    @test_throws ErrorException simpson(Metacommunity(communities,
                                                      GeneralTypes(rand(numspecies,
                                                                        numspecies))))

    @test jaccard([1 0 0; 0 1 1]' / 3)[!, :diversity][1] + 1.0 ≈ 1.0
    @test jaccard(Metacommunity([1 0 1; 0 1 1]' / 4))[!, :diversity] ==
          jaccard([1 0 1; 0 1 1]' / 4)[!, :diversity]
    @test_throws ErrorException jaccard(Metacommunity([1 0 1; 0 1 1]' / 4,
                                                      GeneralTypes(rand(3, 3))))
    @test jaccard([1 0 1; 0 1 1]' / 4)[!, :diversity][1] ≈ 1.0 / 3.0
    @test_throws ErrorException jaccard([1 1 0; 0 1 1; 1 1 1] / 7)
    @test_throws ErrorException jaccard(Metacommunity([1 1 0; 0 1 1; 1 1 1] / 7))

    @test pielou([1, 1])[!, :diversity][1] ≈ 1.0
    @test all(pielou([1 2; 1 2])[!, :diversity] .≈ 1.0)
    @test all(pielou(Metacommunity([1 2; 1 2]))[!, :diversity] .≈ 1.0)

    @test_throws "Can only calculate" gower([1 1 0; 0 1 1; 1 1 1])
    @test gower([1 0; 0 1]).diversity[1] ≈ 1.0
    @test gower([1 1; 1 1]).diversity[1] ≈ 0.0
    @test gower([1 0; 0 0; 0 1; 0 0], countzeros = true).diversity[1] ≈ 0.5
    @test gower([1 0; 0 0; 0 1; 0 0], countzeros = false).diversity[1] ≈ 1.0
end

@testset "Generalised ecological diversities" begin
    @test generalisedrichness(metacommunityDiversity,
                              communities, Z1)[!, :diversity] ≈ [1]
    @test generalisedrichness(subcommunityDiversity,
                              communities)[!, :diversity] ≈
          richness(communities)[!, :diversity]

    @test_throws ErrorException generalisedrichness(individualDiversity,
                                                    communities, Z1)

    @test generalisedshannon(metacommunityDiversity,
                             communities, Z1)[!, :diversity] .+ 1.0 ≈ [1.0]
    @test generalisedshannon(subcommunityDiversity,
                             communities)[!, :diversity] ≈
          shannon(communities)[!, :diversity]
    @test_throws ErrorException generalisedshannon(individualDiversity,
                                                   communities, Z1)

    @test generalisedsimpson(metacommunityDiversity,
                             communities, Z1)[!, :diversity] ≈ [1.0]
    @test generalisedsimpson(subcommunityDiversity,
                             communities)[!, :diversity] ≈
          simpson(communities)[!, :diversity]
    @test_throws ErrorException generalisedsimpson(individualDiversity,
                                                   communities, Z1)

    @test generalisedjaccard([1 0 1; 0 1 1]' / 4, UniqueTypes(3))[!,
                                                                  :diversity][1] ≈
          1.0 / 3.0
    @test generalisedjaccard([1 0 1; 0 1 1]' / 4, ones(Float64, 3, 3))[!,
                                                                       :diversity][1] ≈
          1.0
    @test generalisedjaccard([1 1 1; 1 1 1]' / 6, Matrix(1.0I, 3, 3))[!,
                                                                      :diversity][1] ≈
          1.0

    @test all(generalisedpielou(subcommunityDiversity,
                                [1/6 2/6; 1/6 2/6]).diversity .≈ 1.0)

    mat = reshape(rand(9), 3, 3)
    mat ./= sum(mat)
    mats = Float64[1 0 3; 3 2 0; 4 2 3] ./ 18
    matm = Float64[1 0 1; 1 2 0; 2 2 1] ./ 10
    mat1 = Float64[1 1 1; 1 1 1; 1 1 1] ./ 9
    sim = Float64[1 1 0; 1 1 0; 0 0 1]
    @test all(generalisedpielou(subcommunityDiversity, mat).diversity .≈
              pielou(mat).diversity)
    @test all(generalisedpielou(subcommunityDiversity, mat1).diversity .≈ 1.0)
    @test all(generalisedpielou(metacommunityDiversity, mat1).diversity .≈ 1.0)
    @test all(generalisedpielou(subcommunityDiversity, mat,
                                UniqueTypes(3)).diversity .≈
              pielou(mat).diversity)
    @test_broken all(generalisedpielou(subcommunityDiversity, mats,
                                       sim).diversity .≈ 1.0)
    @test_broken all(generalisedpielou(metacommunityDiversity, matm,
                                       sim).diversity .≈ 1.0)
    @test_throws "Can't calculate Pielou" generalisedpielou(individualDiversity,
                                                            mat, sim)
    @test_throws "function cannot run with" pielou(Metacommunity(mat,
                                                                 GeneralTypes(sim)))
end

@testset "Pairwise dissimilarity" begin
    # 6 types (rows) in 4 subcommunities, and a similarity from points in the
    # plane: 1 - distance / maximum distance
    ab = [3 0 1 2; 0 2 2 0; 1 1 0 4; 0 5 1 0; 2 0 0 1; 1 1 3 0]
    xy = [0.1 0.3; 0.4 0.2; 0.9 0.7; 0.5 0.5; 0.2 0.9; 0.7 0.1]
    dist = [sqrt(sum(abs2, xy[i, :] .- xy[j, :])) for i in 1:6, j in 1:6]
    Z = 1 .- dist ./ maximum(dist)
    lower(D) = [D[i, j] for j in 1:4 for i in (j + 1):4]

    # R adiv::dsimcom(t(ab), Z, type = "dissimilarity"), option "absolute"
    # and "relative"
    adiv = Dict((:sokalsneath, false) =>
                    [0.630204062689044, 0.582711283093493,
                        0.502273143713468, 0.402237249804308,
                        0.629380700688380, 0.744182840732740],
                (:jaccard, false) => [0.460071493514718, 0.411145080139611,
                    0.335356972204401, 0.251750298819423,
                    0.459194395558621, 0.592588527112461],
                (:sorensen, false) => [0.298761592877308, 0.258768170082979,
                    0.201458791227148, 0.144001339539437,
                    0.298022277589717, 0.421048526694661],
                (:ochiai, false) => [0.220053207952419, 0.243928811551166,
                    0.200802053220079, 0.112678666476515,
                    0.232453936301415, 0.413695152327061],
                (:sokalsneath, true) => [0.555769380659253, 0.582711283093493,
                    0.502273143713468, 0.337183100509364,
                    0.564046714055179, 0.744182840732740],
                (:jaccard, true) => [0.384820383404519, 0.411145080139611,
                    0.335356972204401, 0.202778249735525,
                    0.392802968993556, 0.592588527112461],
                (:sorensen, true) => [0.238252377290183, 0.258768170082979,
                    0.201458791227148, 0.112828731182273,
                    0.244402497898829, 0.421048526694661],
                (:ochiai, true) => [0.220053207952419, 0.243928811551166,
                    0.200802053220079, 0.112678666476515,
                    0.232453936301415, 0.413695152327061])
    for relative in (false, true)
        d = pairwisedissimilarity(ab, Z; relative)
        for form in (:jaccard, :sorensen, :ochiai, :sokalsneath)
            @test lower(d[form]) ≈ adiv[(form, relative)] atol = 1e-12
        end
    end

    # R ade4::disc(as.data.frame(ab), as.dist(sqrt(2 * (1 - Z)))), which
    # returns √DISC for Δ = dis², i.e. twice our Δ = 1 - Z
    ade4 = [0.535141765573118, 0.551555371780899, 0.445172500440329,
        0.402206155013214, 0.550920898578258, 0.715389039372716]
    @test sqrt.(2 .* lower(pairwisedissimilarity(ab, Z).rao)) ≈ ade4 atol = 1e-12

    # R vegan::vegdist(t(ab), ...): binary Jaccard and Bray-Curtis, and
    # Morisita-Horn
    d01 = pairwisedissimilarity(ab .> 0)
    @test lower(d01.jaccard) ≈
          [2 / 3, 2 / 3, 0.25, 0.4, 5 / 6, 5 / 6] atol = 1e-12
    @test lower(d01.sorensen) ≈
          [0.5, 0.5, 1 / 7, 0.25, 5 / 7, 5 / 7] atol = 1e-12
    @test lower(pairwisedissimilarity(ab; relative = true).sorensen) ≈
          [0.907827359180688, 0.6, 1 / 3, 0.446964155084126,
        0.843478260869565, 0.888888888888889] atol = 1e-12

    # Rao against its definition, and for distinct types
    p = ab ./ sum(ab, dims = 1)
    Δ = 1 .- Z
    rao = [p[:, x]' * Δ * p[:, y] -
           (p[:, x]' * Δ * p[:, x] + p[:, y]' * Δ * p[:, y]) / 2
           for x in 1:4, y in 1:4]
    @test pairwisedissimilarity(ab, Z).rao ≈ rao atol = 1e-12
    @test pairwisedissimilarity(ab).rao ≈
          [sum(abs2, p[:, x] - p[:, y]) / 2 for x in 1:4, y in 1:4] atol = 1e-12

    # The same answer from proportions, types and a metacommunity; the
    # abundances as given are only up to one common scale, and relative
    # abundances ignore the scale of each subcommunity
    d = pairwisedissimilarity(ab, Z)
    @test pairwisedissimilarity(ab, GeneralTypes(Z)) == d
    @test pairwisedissimilarity(Metacommunity(ab ./ sum(ab),
                                              GeneralTypes(Z))) == d
    @test all(pairwisedissimilarity(10 .* ab, Z)[f] ≈ d[f] for f in keys(d))
    drel = pairwisedissimilarity(ab, Z; relative = true)
    @test all(pairwisedissimilarity(ab .* [1 2 3 4], Z; relative = true)[f] ≈
              drel[f] for f in keys(drel))

    # Symmetric, zero on the diagonal, in [0, 1]; the square root on request
    for D in d
        @test D == D'
        @test all(iszero, diag(D))
        @test all(x -> 0 ≤ x ≤ 1, D)
    end
    dsq = pairwisedissimilarity(ab, Z; squareroot = true)
    @test all(dsq[f] ≈ sqrt.(d[f]) for f in keys(d))

    # A subset of forms, in the order asked for
    @test keys(pairwisedissimilarity(ab, Z; forms = (:rao, :ochiai))) ==
          (:rao, :ochiai)
    @test pairwisedissimilarity(ab, Z; forms = (:ochiai,)).ochiai == d.ochiai

    # An empty subcommunity is NaN, and leaves the other pairs as they were
    ab0 = hcat(ab, zeros(Int, 6))
    d0 = pairwisedissimilarity(ab0, Z)
    @test all(isnan, d0.jaccard[:, 5]) && all(isnan, d0.jaccard[5, :])
    @test d0.jaccard[1:4, 1:4] ≈ d.jaccard

    @test_throws ArgumentError pairwisedissimilarity(ab, Z; forms = (:bray,))
    @test_throws "at least two subcommunities" pairwisedissimilarity(ab[:, 1:1])
    # A similarity that is not positive semi-definite can push s above 1
    notpsd = [1.0 1.0 0.0; 1.0 1.0 1.0; 0.0 1.0 1.0]
    @test_logs (:warn, r"outside \[0, 1\]") match_mode=:any pairwisedissimilarity([1 0;
                                                                                   0 1;
                                                                                   1 0],
                                                                                  notpsd;
                                                                                  forms = (:sorensen,))
end

end
