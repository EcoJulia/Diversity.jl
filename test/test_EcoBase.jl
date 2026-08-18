# SPDX-License-Identifier: BSD-2-Clause

module TestEcoBase
using Test

# Checking EcoBase interface
using Diversity
using EcoBase
# Not in EcoBase's export list, so they must be named explicitly.
using EcoBase: thingkind, thingkindplural, placekind, placekindplural
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

    # views: SpatialEcology's own, on its own assemblage -- asserted rather than just called,
    # which is what it was before.
    va = view(amph, species = 1:10)
    @test counttypes(va) == 10
    @test countsubcommunities(va) == countsubcommunities(amph)

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

@testset "view" begin
    Z = [1.0 0.5 0.0; 0.5 1.0 0.5; 0.0 0.5 1.0]
    types = GeneralTypes(Z, ["ash", "oak", "elm"])
    ab = [0.1 0.2; 0.2 0.1; 0.2 0.2]
    mc = Metacommunity(ab, types, Subcommunities(["north", "south"]))

    # It is a genuine view: a SubArray aliasing the parent, holding a subset that does not sum to
    # one. Only a Metacommunity requires that, and this deliberately is not one.
    v = view(mc, sites = 1:1)
    @test v isa EcoBase.AbstractAssemblage
    @test !(v isa Diversity.AbstractMetacommunity)
    @test occurrences(v) isa SubArray
    @test sum(occurrences(v)) ≈ sum(ab[:, 1])
    @test sum(occurrences(v)) < 1

    # Aliasing, not copying: a change to the parent's abundances shows through. A Metacommunity
    # caches and would not, which is the substantive difference between the two.
    parent = Metacommunity(copy(ab), types, Subcommunities(["north", "south"]))
    alias = view(parent, sites = 1:1)
    before = occurrences(alias)[1, 1]
    getabundance(parent)[1, 1] *= 2
    @test occurrences(alias)[1, 1] ≈ 2 * before

    # The abundances are normalised when read, so the subset measures as a metacommunity in its own
    # right: identical to what a user would have built by hand from the same columns.
    byhand = Metacommunity(ab[:, 1:1] ./ sum(ab[:, 1:1]), types,
                           Subcommunities(["north"]))
    for q in [0, 1, 2, Inf]
        @test meta_gamma(v, q).diversity ≈ meta_gamma(byhand, q).diversity
        @test norm_sub_alpha(v, q).diversity ≈
              norm_sub_alpha(byhand, q).diversity
    end

    # Indices, names and boolean masks all select, because EcoBase's own `asindices` resolves them.
    @test gettypenames(view(mc, species = ["ash", "elm"])) == ["ash", "elm"]
    @test gettypenames(view(mc, species = [1, 3])) == ["ash", "elm"]
    @test gettypenames(view(mc, species = [true, false, true])) ==
          ["ash", "elm"]
    @test getsubcommunitynames(view(mc, sites = ["south"])) == ["south"]

    # Subsetting types materialises the similarity into a GeneralTypes; ash and elm are the pair
    # with no similarity between them, so the submatrix is the identity.
    @test calcsimilarity(gettypes(view(mc, species = [1, 3])), 1) ==
          [1.0 0.0; 0.0 1.0]

    # A dimension that is not restricted keeps its object untouched rather than rebuilding it.
    @test gettypes(view(mc, sites = 1:1)) === gettypes(mc)
    @test getpartition(view(mc, species = 1:2)) === getpartition(mc)
    @test gettypes(view(mc)) === gettypes(mc)

    # UniqueTypes cannot go through the generic default, since its similarity is a UniformScaling.
    umc = Metacommunity(ab)
    @test gettypes(view(umc, species = [1, 3])) isa UniqueTypes
    @test gettypenames(view(umc, species = [1, 3])) == ["1", "3"]

    # An undivided metacommunity stays undivided rather than becoming a Subcommunities of one.
    @test getpartition(view(Metacommunity([0.5, 0.5]), sites = [1])) isa
          Onecommunity

    # Converting gives back a cached Metacommunity measuring the same thing.
    conv = Metacommunity(v)
    @test conv isa Diversity.AbstractMetacommunity
    @test meta_gamma(conv, 1).diversity ≈ meta_gamma(v, 1).diversity

    # It prints as one of ours, with the singular where it belongs -- EcoBase's own `show` for a
    # generic assemblage would say "1 subcommunities" here.
    out = sprint(show, v)
    @test occursin("with 3 species in 1 subcommunity ", out)

    # A subset names its units from its own types and partition rather than falling back to
    # EcoBase's "thing" and "place". Both plurals need more than one of something to appear, which
    # is why they are asserted here as well as through `show`.
    whole = view(mc)
    @test thingkind(whole) == "species"
    @test thingkindplural(whole) == "species"
    @test placekind(whole) == "subcommunity"
    @test placekindplural(whole) == "subcommunities"
    @test occursin("with 3 species in 2 subcommunities", sprint(show, whole))

    # This is what the missing `view` was blocking in EcoBase itself.
    @test EcoBase.cooccurring(mc, [1, 2]) == [true, true]
end

@testset "view unlocks SpatialEcology grouping" begin
    # `groupsites` and `groupspecies` are typed on EcoBase.AbstractAssemblage rather than on
    # SpatialEcology's own types, so implementing `view` is all that was needed to make them work
    # on a metacommunity. They were a MethodError before.
    mc = Metacommunity([0.1 0.2 0.1; 0.2 0.1 0.1; 0.1 0.05 0.05],
                       UniqueTypes(["a", "b", "c"]),
                       Subcommunities(["x", "y", "z"]))
    groups = groupsites(mc, ["left", "left", "right"])
    @test length(groups) == 2
    @test all(g -> g isa EcoBase.AbstractAssemblage, groups)
    @test getsubcommunitynames.(groups) == [["x", "y"], ["z"]]
    # Each group is measurable, and the two-site group is the one with two subcommunities.
    @test nrow(norm_sub_alpha(groups[1], 1)) == 2
    @test nrow(norm_sub_alpha(groups[2], 1)) == 1

    species = groupspecies(mc, ["plant", "plant", "animal"])
    @test length(species) == 2
    @test sort(vcat(gettypenames.(species)...)) == ["a", "b", "c"]
end

end
