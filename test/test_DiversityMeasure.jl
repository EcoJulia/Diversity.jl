module TestDiversityMeasure
using Compat.Test
using Compat.LinearAlgebra
using Compat

using Diversity
using Diversity.ShortNames
using DataFrames

pop = [3, 3, 4]
pop = pop / sum(pop)
oc = Onecommunity()
sim = [1.0 0 0; 1.0 1.0 0.0; 1.0 1.0 1.0]
ms = GeneralTypes(sim)
meta = Metacommunity(pop, ms, oc)
meta1 = Metacommunity(pop)

ab3 = [1.0 2; 3 0; 0 4]'
ab3 = ab3 / sum(ab3)
sc = Subcommunities(3)
sp = Species(2)
meta2 = Metacommunity(ab3, sp, sc)
nab = NormalisedAlpha(meta2)

@testset "Diversity measures" begin
    diversities = [RawAlpha, NormalisedAlpha, RawBeta, NormalisedBeta,
                   RawRho, NormalisedRho, Gamma]
    shortds = [α, ᾱ, β, β̄, ρ, ρ̄, Γ]
    chars = ["α", "ᾱ", "β", "β̄", "ρ", "ρ̄", "γ"]
    asciis = ["RawAlpha", "NormalisedAlpha",
              "RawBeta", "NormalisedBeta",
              "RawRho", "NormalisedRho", "Gamma"]
    fulls = ["raw alpha diversity", "normalised alpha diversity",
             "distinctiveness", "effective number of subcommunities",
             "redundancy", "representativeness", "gamma diversity"]
    for i in 1:length(diversities)
        @test diversities[i] == shortds[i]
        @test getName(diversities[i](meta)) == chars[i]
        @test getASCIIName(diversities[i](meta2)) == asciis[i]
        @test getFullName(diversities[i](meta1)) == fulls[i]
    end
end

numbers = [1., 2, 4, 8, 16];
numspecies = 100;
fragments = rand(numspecies);
weights = rand(numspecies);
weights /= sum(weights);
Z1 = ones(typeof(weights[1]), (length(weights), length(weights)));
numcommunities = 8;
manyweights = rand(numspecies, numcommunities);
manyweights *= Diagonal(reshape(mapslices(v -> 1. / sum(v), manyweights;
                                          dims=1), numcommunities));

@testset "inddiv / subdiv / metadiv" begin
    @test individualDiversity(nab, 0)[:diversity] ≈ inddiv(nab, 0)[:diversity]
    @test individualDiversity(nab)(1)[:diversity] ≈ inddiv(nab, 1)[:diversity]
    idnab = inddiv(nab, [2, 3])
    @test idnab[isapprox.(collect(idnab[:q]), 2), :diversity] ≈ inddiv(nab, 2)[:diversity]
    allid = inddiv(meta2, Inf)
    @test allid[allid[:measure] .== "RawAlpha", :diversity] ≈ inddiv(RawAlpha(meta2), Inf)[:diversity]

    @test subcommunityDiversity(nab, Inf)[:diversity] ≈ subdiv(nab, Inf)[:diversity]
    sdnab = subdiv(nab, [4, 5])
    @test sdnab[isapprox.(collect(sdnab[:q]), 4), :diversity] ≈ subdiv(nab, 4)[:diversity]
    allsd = subdiv(meta2, Inf)
    @test allsd[allsd[:measure] .== "NormalisedAlpha", :diversity] ≈ subdiv(nab, Inf)[:diversity]
    allmd = subdiv(meta1, Inf)
    @test allmd[allmd[:measure] .== "NormalisedAlpha", :diversity] ≈ metadiv(ᾱ(meta1), Inf)[:diversity]

    scg = subcommunityDiversity(Gamma(meta2))
    @test scg(1)[:diversity] ≈ scg(1.0)[:diversity]

    communities = rand(numspecies, numcommunities);
    communities /= sum(communities);
    @test subdiv(NormalisedAlpha(Metacommunity(communities)), 0)[:diversity] ≈ numspecies * ones(size(communities, 2))
    @test subdiv(NormalisedAlpha(Metacommunity(communities)), [0])[:diversity] ≈ numspecies * ones(size(communities, 2))
    qs = [0, 1, 2, Inf]
    sna = subdiv(NormalisedAlpha(Metacommunity(communities, Z1)), qs)
    @test nrow(sna) == length(qs) * size(communities, 2)
    for q in qs
        @test sna[isapprox.(sna[:q], q), :diversity] ≈ ones(size(communities, 2))
    end
    @test subdiv(RawAlpha(Metacommunity(communities)), 0)[:diversity] ≈ numspecies * vec(mapslices(v -> 1. / sum(v), communities; dims=1))

    even = ones((numspecies, numcommunities)) / (numspecies * numcommunities);
    qs = [0, 1, 2, 3, 4, 5, 6, Inf];
    @test metadiv(NormalisedAlpha(Metacommunity(even)), qs)[:diversity] ≈ numspecies * ones(length(qs))
    @test metadiv(RawAlpha(Metacommunity(even)), qs)[:diversity] ≈ numspecies * numcommunities * ones(length(qs))
    md2 = metadiv(meta2, Inf)
    @test md2[md2[:measure] .== "Gamma", :diversity] ≈ metadiv(Gamma(meta2), Inf)[:diversity]

    probs = reshape(mapslices(sum, communities; dims=2),
                    size(communities, 1));
    @test metadiv(Gamma(Metacommunity(communities)), qs)[:diversity] ≈ qD(probs, qs)
    @test metadiv(Gamma(Metacommunity(communities, Z1)), qs)[:diversity] ≈ qDZ(probs, qs, Z1)

    Z = rand(numspecies, numspecies)
    @test metadiv(Gamma(Metacommunity(communities, Z)), qs)[:diversity] ≈ qDZ(probs, qs, Z)

    colweights = rand(numcommunities)
    colweights /= sum(colweights)
    allthesame = probs * colweights'
    @test metadiv(RawBeta(Metacommunity(allthesame, Z)), qs)[:diversity] ≈ 1.0 ./ qD(colweights, 2 .- qs)
    @test metadiv(NormalisedBeta(Metacommunity(allthesame, Z)), qs)[:diversity] ≈ ones(length(qs))
    @test metadiv(NormalisedRho(Metacommunity(allthesame, Z)), qs)[:diversity] ≈ ones(length(qs))
    @test metadiv(RawRho(Metacommunity(allthesame, Z)), qs)[:diversity] ≈ qD(colweights, qs)

    communitylist = rand(1:numcommunities, numspecies)
    distinct = zeros(Float64, (numspecies, numcommunities))
    for i in 1:numspecies
        distinct[i, communitylist[i]] = weights[i]
    end

    @test metadiv(RawRho(Metacommunity(distinct)), qs)[:diversity] ≈ ones(length(qs))
    subnr = subdiv(NormalisedRho(Metacommunity(distinct)), qs)
    for q in qs
        @test subnr[isapprox.(subnr[:q], q), :diversity] ≈
              vec(Compat.sum(distinct; dims=1))
    end
    @test metadiv(NormalisedBeta(Metacommunity(distinct)), qs)[:diversity] ≈
        qD(reshape(Compat.sum(distinct; dims=1), numcommunities), qs)
    @test metadiv(RawBeta(Metacommunity(distinct)), qs)[:diversity] ≈ ones(length(qs))

    # many (unexported!) diversity levels not yet implemented
    @test_throws ErrorException Diversity.communityDiversity(nab)
end

end
