module TestDiversityMeasure
using Diversity

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

oc = Onecommunity([3, 3, 4])
sim = [1.0 0 0; 1.0 1.0 0.0; 1.0 1.0 1.0]
ms = MatrixSimilarity(sim)
sup = Supercommunity(oc, ms)
sup3 = Supercommunity(oc)

ab3 = [1 2; 3 0; 0 4]'
sc = Subcommunities(ab3)
sp = Species()
sup2 = Supercommunity(sc, sp)
nab = NormalisedAlpha(sup2)
    

@testset "Diversity measures" begin
    nr = NormalisedRho(sup)
    @test getName(nr) == "ρ̄"
    @test getASCIIName(nr) == "rho bar"
    @test getFullName(nr) == "representativeness"
    
    @test getName(nab) == "ᾱ"
    @test getASCIIName(nab) == "alpha bar"
    @test getFullName(nab) == "normalised alpha diversity"
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


# Sub-community alpha diversities

@testset "inddiv / subdiv / superdiv" begin
    @test individualDiversity(nab, 0) ≈ inddiv(nab, 0)
    @test individualDiversity(nab)(1) ≈ inddiv(nab, 1)
    @test subcommunityDiversity(nab, Inf) ≈ subdiv(nab, Inf)

    scg = subcommunityDiversity(Gamma(sup2))
    @test scg(1) ≈ scg(1.0)
    
    communities = rand(numspecies, numcommunities);
    communities /= sum(communities);
    @test subdiv(NormalisedAlpha(Supercommunity(communities)), 0) ≈ numspecies * ones(size(communities, 2))
    @test subdiv(NormalisedAlpha(Supercommunity(communities)),
                 [0])[1] ≈ numspecies * ones(size(communities, 2))
    sna = subdiv(NormalisedAlpha(Supercommunity(communities, Z1)), [0, 1, 2, Inf])
    @test length(sna) == 4
    for i in eachindex(sna)
        @test sna[i] ≈ ones(size(communities, 2))
    end
    
    @test subdiv(RawAlpha(Supercommunity(communities)), 0) ≈ numspecies * vec(mapslices(v -> 1. / sum(v), communities, 1))
    
    even = ones((numspecies, numcommunities)) / (numspecies * numcommunities);
    qs = [0, 1, 2, 3, 4, 5, 6, Inf];
    @test superdiv(NormalisedAlpha(Supercommunity(even)), qs) ≈ numspecies * ones(length(qs))
    @test superdiv(RawAlpha(Supercommunity(even)), qs) ≈ numspecies * numcommunities * ones(length(qs))
    
    probs = reshape(mapslices(sum, communities, 2), (size(communities, 1)));
    @test superdiv(Gamma(Supercommunity(communities)), qs) ≈ qD(probs, qs)
    @test superdiv(Gamma(Supercommunity(communities, Z1)), qs) ≈ qDZ(probs, qs, Z1)
    
    Z = rand(numspecies, numspecies);
    @test superdiv(Gamma(Supercommunity(communities, Z)), qs) ≈ qDZ(probs, qs, Z)
    
    colweights = rand(numcommunities);
    colweights /= sum(colweights);
    allthesame = probs * colweights';
    @test superdiv(RawBeta(Supercommunity(allthesame, Z)), qs) ≈ 1.0 ./ qD(colweights, 2 - qs)
    @test superdiv(NormalisedBeta(Supercommunity(allthesame, Z)), qs) ≈ ones(length(qs))
    @test superdiv(NormalisedRho(Supercommunity(allthesame, Z)), qs) ≈ ones(length(qs))
    @test superdiv(RawRho(Supercommunity(allthesame, Z)), qs) ≈ qD(colweights, qs)
    
    communitylist = rand(1:numcommunities, numspecies)
    distinct = zeros(Float64, (numspecies, numcommunities))
    for i in 1:numspecies
        distinct[i, communitylist[i]] = weights[i]
    end
    
    @test superdiv(RawRho(Supercommunity(distinct)), qs) ≈ ones(length(qs))
    subnr = subdiv(NormalisedRho(Supercommunity(distinct)), qs)
    for i in eachindex(qs)
        @test subnr[i] ≈ vec(sum(distinct, 1))
    end
    @test superdiv(NormalisedBeta(Supercommunity(distinct)), qs) ≈ qD(reshape(sum(distinct, 1), numcommunities), qs)
    @test superdiv(RawBeta(Supercommunity(distinct)), qs) ≈ ones(length(qs))
end
    
end
