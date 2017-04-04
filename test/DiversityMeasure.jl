module TestDiversityMeasure
using Diversity
using Diversity.ShortNames
using Base.Test

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
manyweights *= diagm(reshape(mapslices(v -> 1. / sum(v), manyweights, 1),
                             (numcommunities)));

@testset "inddiv / subdiv / metadiv" begin
    @test individualDiversity(nab, 0) ≈ inddiv(nab, 0)
    @test individualDiversity(nab)(1) ≈ inddiv(nab, 1)
    @test inddiv(nab, [2, 3])[1] ≈ inddiv(nab, 2)
    @test inddiv(meta2, Inf)[1] ≈ inddiv(RawAlpha(meta2), Inf)

    @test subcommunityDiversity(nab, Inf) ≈ subdiv(nab, Inf)
    @test subdiv(nab, [4, 5])[1] ≈ subdiv(nab, 4)
    @test subdiv(meta2, Inf)[2] ≈ subdiv(nab, Inf)
    @test metadiv(meta1, Inf)[2] ≈ metadiv(ᾱ(meta1), Inf)

    scg = subcommunityDiversity(Gamma(meta2))
    @test scg(1) ≈ scg(1.0)

    communities = rand(numspecies, numcommunities);
    communities /= sum(communities);
    @test subdiv(NormalisedAlpha(Metacommunity(communities)), 0) ≈ numspecies * ones(size(communities, 2))
    @test subdiv(NormalisedAlpha(Metacommunity(communities)),
                 [0])[1] ≈ numspecies * ones(size(communities, 2))
    sna = subdiv(NormalisedAlpha(Metacommunity(communities, Z1)), [0, 1, 2, Inf])
    @test length(sna) == 4
    for i in eachindex(sna)
        @test sna[i] ≈ ones(size(communities, 2))
    end 
    @test subdiv(RawAlpha(Metacommunity(communities)), 0) ≈ numspecies * vec(mapslices(v -> 1. / sum(v), communities, 1))

    even = ones((numspecies, numcommunities)) / (numspecies * numcommunities);
    qs = [0, 1, 2, 3, 4, 5, 6, Inf];
    @test metadiv(NormalisedAlpha(Metacommunity(even)), qs) ≈ numspecies * ones(length(qs))
    @test metadiv(RawAlpha(Metacommunity(even)), qs) ≈ numspecies * numcommunities * ones(length(qs))
    @test metadiv(meta2, Inf)[7] ≈ metadiv(Gamma(meta2), Inf)

    probs = reshape(mapslices(sum, communities, 2), (size(communities, 1)));
    @test metadiv(Gamma(Metacommunity(communities)), qs) ≈ qD(probs, qs)
    @test metadiv(Gamma(Metacommunity(communities, Z1)), qs) ≈ qDZ(probs, qs, Z1)

    Z = rand(numspecies, numspecies);
    @test metadiv(Gamma(Metacommunity(communities, Z)), qs) ≈ qDZ(probs, qs, Z)

    colweights = rand(numcommunities);
    colweights /= sum(colweights);
    allthesame = probs * colweights';
    @test metadiv(RawBeta(Metacommunity(allthesame, Z)), qs) ≈ 1.0 ./ qD(colweights, 2 - qs)
    @test metadiv(NormalisedBeta(Metacommunity(allthesame, Z)), qs) ≈ ones(length(qs))
    @test metadiv(NormalisedRho(Metacommunity(allthesame, Z)), qs) ≈ ones(length(qs))
    @test metadiv(RawRho(Metacommunity(allthesame, Z)), qs) ≈ qD(colweights, qs)

    communitylist = rand(1:numcommunities, numspecies)
    distinct = zeros(Float64, (numspecies, numcommunities))
    for i in 1:numspecies
        distinct[i, communitylist[i]] = weights[i]
    end

    @test metadiv(RawRho(Metacommunity(distinct)), qs) ≈ ones(length(qs))
    subnr = subdiv(NormalisedRho(Metacommunity(distinct)), qs)
    for i in eachindex(qs)
        @test subnr[i] ≈ vec(sum(distinct, 1))
    end
    @test metadiv(NormalisedBeta(Metacommunity(distinct)), qs) ≈ qD(reshape(sum(distinct, 1), numcommunities), qs)
    @test metadiv(RawBeta(Metacommunity(distinct)), qs) ≈ ones(length(qs))

    # many (unexported!) diversity levels not yet implemented
    @test_throws ErrorException Diversity.communityDiversity(nab)
end

end
