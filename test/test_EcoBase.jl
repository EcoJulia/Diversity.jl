# SPDX-License-Identifier: BSD-2-Clause

module TestEcoBase
using Test

# Checking EcoBase interface
using Diversity
using EcoBase
using SpatialEcology
using CSV
using DataFrames
using Plots

numspecies = 10
numcommunities = 8
manyweights = rand(numspecies, numcommunities)
manyweights /= sum(manyweights)

# Subtypes that implement neither the Diversity API nor the EcoBase interface — the case the guards
# in src/EcoBase.jl exist for. Fully qualified because SpatialEcology is loaded here too.
struct BarePartition <: Diversity.AbstractPartition{Nothing} end
struct BareTypes <: Diversity.AbstractTypes end
struct BareMC <:
       Diversity.AbstractMetacommunity{Float64, Matrix{Float64},
                                       Matrix{Float64}, BareTypes,
                                       BarePartition} end

# Two types that name themselves but supply no similarity matrix, differing only in whether they
# claim to have similarity at all — which is what decides whether the identity default is right.
struct ClaimsSimilarity <: Diversity.AbstractTypes end
struct DeclaresNoSimilarity <: Diversity.AbstractTypes end
Diversity.API._gettypenames(::ClaimsSimilarity, ::Bool) = ["a", "b"]
Diversity.API._gettypenames(::DeclaresNoSimilarity, ::Bool) = ["a", "b"]
Diversity.API._hassimilarity(::DeclaresNoSimilarity) = false

@testset "EcoBase interface" begin
    species = map(n -> "Species $n", 1:numspecies)
    communities = map(n -> "SC $n", 1:numcommunities)
    ut = UniqueTypes(species)
    sc = Subcommunities(communities)
    mc = Metacommunity(manyweights, ut, sc)

    @test all(thingnames(mc) .== species)
    @test all(placenames(mc) .== communities)
    @test all(occurrences(mc) .≈ manyweights)
    @test all(richness(mc) .== repeat([numspecies], inner = numcommunities))
    @test all(occupancy(mc) .== repeat([numcommunities], inner = numspecies))
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
    amphdat = CSV.read(joinpath(dirname(pathof(SpatialEcology)), "..", "data",
                                "amph_Europe.csv"), DataFrame)
    amph = Assemblage(amphdat[!, 4:end], amphdat[!, 1:3], sitecolumns = false)

    @test typeof(amph) ==
          Assemblage{Bool, SpatialEcology.Locations{SpatialEcology.GridData}}

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
    amp2 = coarsen(amph, 2)
    @test sum(richness(amp2)) == 2862
    @test nsites(amp2) == 285

    pointamph = Assemblage(amphdat[!, 4:end], amphdat[!, 1:3],
                           sitecolumns = false,
                           cdtype = SpatialEcology.pointdata)
    amp3 = coarsen(pointamph, 2)
    @test richness(amp3) == richness(amp2)
    rich2 = metadiv(Gamma(amp2), 0)
    rich3 = metadiv(Gamma(amp3), 0)
    @test rich2.diversity[1] == rich3.diversity[1] ==
          meta_gamma(Metacommunity(amp3), 0).diversity[1]
    @test getaddedoutput(amph) === nothing
    @test getaddedoutput(gettypes(amph)) === nothing

    # This is the route to the spatial diversity maps the framework was designed to produce.
    @test plot(norm_sub_rho(amph, 1), amph) isa Plots.Plot
    @test plot(sub_gamma(amph, 0), amph) isa Plots.Plot
    @test nrow(norm_sub_rho(amph, 1)) == countsubcommunities(amph)
end

@testset "EcoBase bridge does not recurse" begin
    # A subtype that implements neither the Diversity API nor the EcoBase interface must say so,
    # not overflow the stack. Regression for the two-way bridge in src/EcoBase.jl.
    @test_throws ErrorException Diversity.API._getpartition(BareMC())
    @test_throws ErrorException Diversity.API._gettypes(BareMC())
    @test_throws ErrorException Diversity.API._getabundance(BareMC(), true)
    @test_throws ErrorException Diversity.API._getsubcommunitynames(BarePartition())
    @test_throws ErrorException Diversity.API._gettypenames(BareTypes(), true)

    # …and the reverse bridge must still work for a genuine non-Diversity assemblage. This overlaps
    # the SpatialEcology testset above deliberately: that one is the real coverage, on real data,
    # while this hand-built four-site case sits beside the guards and shows what they must not break.
    asm = Assemblage([1 0 1 1; 0 1 1 0; 1 1 0 1],
                     Float64[1 1; 2 1; 1 2; 2 2],
                     ["s1", "s2", "s3", "s4"], ["a", "b", "c"])
    @test size(norm_sub_alpha(asm, 1.0)) == (4, 8)
    @test size(meta_gamma(asm, 1.0)) == (1, 8)
end

@testset "Similarity claimed but not implemented" begin
    # `_calcsimilarity` fails the other way: the EcoBase fallback answers it with an identity
    # matrix rather than recursing, so a type that claims similarity and never supplied one used
    # to get silently wrong numbers. It is the claim that is refused, not the omission —
    # declaring no similarity still earns the identity default.
    @test_throws ErrorException calcsimilarity(ClaimsSimilarity(), 1)
    @test calcsimilarity(DeclaresNoSimilarity(), 1) == [1.0 0.0; 0.0 1.0]

    # And the types the package ships are untouched, similarity or not.
    @test calcsimilarity(UniqueTypes(2), 1) == [1.0 0.0; 0.0 1.0]
    Z = [1.0 0.3; 0.3 1.0]
    @test calcsimilarity(GeneralTypes(Z), 1) == Z
end

end
