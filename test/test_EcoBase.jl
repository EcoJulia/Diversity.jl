module TestEcoBase
using Test

# Checking EcoBase interface
using Diversity
using EcoBase
using SpatialEcology
using CSV
using DataFrames

numspecies = 10;
numcommunities = 8;
manyweights = rand(numspecies, numcommunities);
manyweights /= sum(manyweights);

@testset "EcoBase interface" begin
    species = map(n -> "Species $n", 1:numspecies)
    communities = map(n -> "SC $n", 1:numcommunities)
    ut = UniqueTypes(species)
    sc = Subcommunities(communities)
    mc = Metacommunity(manyweights, ut, sc)

    @test all(thingnames(mc) .== species)
    @test all(placenames(mc) .== communities)
    @test all(occurrences(mc) .≈ manyweights)
    @test all(richness(mc) .== repeat([numspecies], inner=numcommunities))
    @test all(occupancy(mc) .== repeat([numcommunities], inner=numspecies))
    fewerweights = deepcopy(manyweights)
    fewerweights[1, 1] = 0
    fewerweights /= sum(fewerweights)
    fmc = Metacommunity(fewerweights, ut, sc)

    @test noccupied(fmc) == numcommunities
    @test noccurring(fmc) == numspecies
    @test noccupied(fmc, 1) == numcommunities - 1
    @test noccurring(fmc, 1) == numspecies - 1
    @test nthings(fmc) == numspecies
    @test nplaces(fmc) == numcommunities
    coord = coordinates(getpartition(Metacommunity(reshape(1:81, 9, 9))))
    @test coord[1, 1] ≈ 1.0
    @test coord[9, 2] ≈ 3.0
end

@testset "SpatialEcology.Assemblage" begin
    amphdat = CSV.read(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv"), DataFrame)
    amph = Assemblage(amphdat[!,4:end], amphdat[!,1:3], sitecolumns = false)

    @test typeof(amph) == Assemblage{Bool,SpatialEcology.Locations{SpatialEcology.GridData}}

    # accesseors
    @test extrema(richness(amph)) == (1, 20)
    @test countsubcommunities(amph) == 1010
    @test length(getsubcommunitynames(amph)) == 1010
    @test counttypes(amph) == 73
    @test length(gettypenames(amph)) == 73
    @test occurring(amph, 718) == [15]
    @test occurring(amph, 718:729) == [15, 46, 53, 56]
    @test occupied(amph, "Pleurodeles_waltl")[2] == 14
    @test occupied(amph, ["Pleurodeles_waltl", "Salamandra_corsica"])[50] == 885
    @test occupancy(amph)[1] == 353

    # views
    va = view(amph, species = 1:10)

    #operations
    amp2 = aggregate(amph, 2)
    @test sum(richness(amp2)) == 2862
    @test nsites(amp2) == 285

    pointamph = Assemblage(amphdat[!,4:end], amphdat[!,1:3], sitecolumns = false, cdtype = SpatialEcology.pointdata)
    amp3 = aggregate(pointamph, 2)
    @test richness(amp3) == richness(amp2)
    rich2 = metadiv(Gamma(amp2), 0)
    rich3 = metadiv(Gamma(amp3), 0)
    @test rich2.diversity[1] == rich3.diversity[1] == meta_gamma(Metacommunity(amp3), 0).diversity[1]
    @test getaddedoutput(amph) === nothing
    @test getaddedoutput(gettypes(amph)) === nothing
end

end
