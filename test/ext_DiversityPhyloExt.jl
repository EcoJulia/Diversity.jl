# SPDX-License-Identifier: BSD-2-Clause

module TestPhylogenetics
using Test

using Phylo
using Diversity
using Diversity.Ecology: faith_pd, generalisedfaith_pd
using EcoBase: thingkind, thingkindplural, placekind

# A phylogenetic type that is *not* this extension's `PhyloBranches` — enough of the API to build a
# metacommunity from, and no more. Faith's PD must decline to run on it.
struct OtherPhyloTypes <: Diversity.AbstractPhyloTypes{Nothing} end
Diversity.API._gettypenames(::OtherPhyloTypes, ::Bool) = ["a", "b"]
Diversity.API._calcsimilarity(::OtherPhyloTypes, ::Real) = [1.0 0.5; 0.5 1.0]

@testset "Phylo" begin
    species = ["Dog", "Human", "Cat"]
    abund = [0.4, 0.3, 0.3]
    nt = RootedTree(species)
    n = createnode!(nt)
    createbranch!(nt, n, species[1], 1.0)
    createbranch!(nt, n, species[2], 1.0)
    r = createnode!(nt)
    createbranch!(nt, r, n, 1.0)
    createbranch!(nt, r, species[3], 2.0)
    ph = PhyloBranches(nt)
    leafnames = gettypenames(ph, true)
    @test species == gettypenames(ph, true)
    metaphylo = Metacommunity(abund, ph)
    @test gettypenames(metaphylo, true) == species
    @test getabundance(metaphylo, true) ≈ abund
    @test getabundance(metaphylo) ≈ [0.2, 0.2, 0.15, 0.15, 0.3]
    @test getordinariness!(metaphylo) ≈ [0.4, 0.7, 0.3, 0.7, 0.3]
    @test calcsimilarity(ph, metaphylo.scale) * metaphylo.processedabundances ≈
          getordinariness!(metaphylo)
    @test meta_gamma(metaphylo, 0)[!, :diversity] == [2.5]
    @test sub_gamma(metaphylo, 0)[!, :diversity] == [2.5]

    # Iterating a metacommunity whose processed types (branches) outnumber its raw ones (species).
    # `test_Iterators.jl` covers the iterators themselves without needing `Phylo`; this is the case
    # only a phylogeny produces, so it lives here where `Phylo` is loaded anyway.
    ti = TypeIterator(metaphylo)
    si = SubcommunityIterator(metaphylo)
    @test length(ti) == counttypes(metaphylo, false) == 5
    @test length(ti) > counttypes(metaphylo, true)
    @test length(si) == countsubcommunities(metaphylo)
    @test all(reduce(+, ti) .≈ getweight(metaphylo))
    @test sum(sum, ti) ≈ 1.0

    tsph = PhyloBranches(TreeSet(Dict("tree" => nt)))
    @test species == gettypenames(tsph, true)
    tsmetaphylo = Metacommunity(abund, tsph)
    @test species == gettypenames(tsmetaphylo, true)
    @test meta_gamma(tsmetaphylo, 0).treename == ["tree"]
    @test subdiv(Gamma(tsmetaphylo), 0).treename == ["tree"]
    @test metadiv(Gamma(tsmetaphylo), 0).treename == ["tree"]
    @test all(inddiv(Gamma(tsmetaphylo), 0).treename .== "tree")

    # Note: The units of a phylogenetic metacommunity are *branches*, and this is how a reader is told
    # — `PhyloBranches` is opinionated about that, so it answers EcoBase's naming hooks itself.
    @test thingkind(metaphylo) == "branch"
    @test thingkindplural(metaphylo) == "branches"   # not EcoBase's default "branchs"
    @test placekind(metaphylo) == "subcommunity"
    out = sprint(show, metaphylo)
    @test occursin("with 5 branches in 1 subcommunity", out)
    @test occursin("Branch names:", out)             # singular in the heading

    # Translating to a GeneralTypes metacommunity has to carry the *scaled* Zmatrix and the
    # *branch* abundances together, or the phylogeny's numbers do not survive — which is the case
    # that makes the scale argument to `_calcsimilarity` load-bearing here and nowhere else.
    conv = Metacommunity(metaphylo)
    @test gettypes(conv) isa GeneralTypes
    @test gettypenames(conv) == gettypenames(metaphylo)
    for q in [0, 1, 2, Inf]
        @test norm_sub_alpha(conv, q).diversity ≈
              norm_sub_alpha(metaphylo, q).diversity
        @test norm_sub_rho(conv, q).diversity ≈
              norm_sub_rho(metaphylo, q).diversity
        @test meta_gamma(conv, q).diversity ≈ meta_gamma(metaphylo, q).diversity
    end
end

@testset "Faith's PD" begin
    # Faith's PD is the total length of the branches spanned by the types present, with no
    # normalisation — so the numbers here are read straight off the tree, not off the framework.
    species = ["Dog", "Human", "Cat"]
    nt = RootedTree(species)              # 1 + 1 + 1 + 2 = 5.0 of branch, ultrametric
    n = createnode!(nt)
    createbranch!(nt, n, species[1], 1.0)
    createbranch!(nt, n, species[2], 1.0)
    r = createnode!(nt)
    createbranch!(nt, r, n, 1.0)
    createbranch!(nt, r, species[3], 2.0)
    ph = PhyloBranches(nt)

    @test generalisedfaith_pd(metacommunityDiversity,
                              Metacommunity([0.4, 0.3, 0.3],
                                            ph)).diversity[1] ≈ 5.0

    @test generalisedfaith_pd(metacommunityDiversity,
                              Metacommunity([0.1, 0.1, 0.8],
                                            ph)).diversity[1] ≈ 5.0
    # Drop Dog and the branch to it goes too, but the one it shared with Human stays: 5 - 1 = 4.
    @test generalisedfaith_pd(metacommunityDiversity,
                              Metacommunity([0.0, 0.5, 0.5],
                                            ph)).diversity[1] ≈ 4.0

    # Non-ultrametric, so the scale is a genuinely weighted mean: 1 + 3 + 1 + 2 = 7.0
    nt2 = RootedTree(species)
    n2 = createnode!(nt2)
    createbranch!(nt2, n2, species[1], 1.0)
    createbranch!(nt2, n2, species[2], 3.0)
    r2 = createnode!(nt2)
    createbranch!(nt2, r2, n2, 1.0)
    createbranch!(nt2, r2, species[3], 2.0)
    @test generalisedfaith_pd(metacommunityDiversity,
                              Metacommunity([0.4, 0.3, 0.3],
                                            PhyloBranches(nt2))).diversity[1] ≈
          7.0

    # Per subcommunity, each in isolation: Dog+Human spans 1 + 1 + 1, Human+Cat spans 1 + 1 + 2.
    mc = Metacommunity([0.4 0.0; 0.1 0.2; 0.0 0.3], ph)
    @test faith_pd(mc).diversity ≈ [3.0, 4.0]
    @test faith_pd(mc).diversity ==
          generalisedfaith_pd(subcommunityDiversity, mc).diversity

    # There is no q to ask for, so the column is dropped rather than given a value.
    pd = faith_pd(mc)
    @test "q" ∉ names(pd)
    @test all(pd.measure .== "Faith's PD")
    @test_throws ErrorException generalisedfaith_pd(individualDiversity, mc)

    # It is only defined for phylogenetic types — there is no PD without a tree.
    @test_throws MethodError faith_pd(Metacommunity([0.5, 0.5]))

    # Narrower than that, in fact: only for the `PhyloBranches` this extension supplies, not for
    # `AbstractPhyloTypes` at large. The scale is a total branch length only because that type's
    # `_calcabundance` makes it one; another phylogenetic type may process abundances differently,
    # and would get a confidently wrong number rather than a refusal if the signature were widened.
    # Note: Reached through `get_extension` on purpose: the bare name `PhyloBranches` here is the
    # *abstract* one the parent exports, not the concrete struct the signature is written against.
    concrete = Base.get_extension(Diversity, :DiversityPhyloExt).PhyloBranches
    @test concrete <: Diversity.PhyloBranches <: Diversity.AbstractPhyloTypes
    @test OtherPhyloTypes <: Diversity.AbstractPhyloTypes
    @test !(OtherPhyloTypes <: concrete)
    @test_throws MethodError faith_pd(Metacommunity([0.5, 0.5],
                                                    OtherPhyloTypes()))
end

@testset "view over branches" begin
    species = ["Dog", "Human", "Cat"]
    nt = RootedTree(species)
    n = createnode!(nt)
    createbranch!(nt, n, species[1], 1.0)
    createbranch!(nt, n, species[2], 1.0)
    r = createnode!(nt)
    createbranch!(nt, r, n, 1.0)
    createbranch!(nt, r, species[3], 2.0)
    ph = PhyloBranches(nt)
    mc = Metacommunity([0.4 0.2; 0.1 0.1; 0.1 0.1], ph)

    # `species` selects things, and for this type a thing is a branch -- there are five of them for
    # three species. That is the whole point of the type, and the keyword name is EcoBase's.
    @test counttypes(mc) == 5
    @test length(gettypenames(mc)) == 5
    @test gettypenames(view(mc, species = [1, 2, 3])) == gettypenames(mc)[1:3]

    # Restricting only the subcommunities leaves the phylogeny alone, so nothing is lost.
    sites = view(mc, sites = 1:1)
    @test gettypes(sites) === gettypes(mc)
    @test getdiversityname(gettypes(sites)) == "Phylogenetic Branch"
    @test thingkind(sites) == "branch"

    # Restricting branches cannot: an arbitrary subset of branches is not a tree, so it becomes a
    # GeneralTypes carrying the *scaled* similarity submatrix. That is what keeps the numbers
    # right, and it is why the result reports itself as arbitrary rather than phylogenetic.
    sub = view(mc, species = [1, 2, 3])
    @test gettypes(sub) isa GeneralTypes
    @test getdiversityname(gettypes(sub)) == "Arbitrary Z"

    scaled = calcsimilarity(ph, Diversity.API._getscale(mc))
    ab = getabundance(mc)[1:3, :]
    byhand = Metacommunity(ab ./ sum(ab),
                           GeneralTypes(scaled[1:3, 1:3],
                                        gettypenames(mc)[1:3]))
    for q in [0, 1, 2, Inf]
        @test meta_gamma(sub, q).diversity ≈ meta_gamma(byhand, q).diversity
        @test norm_sub_alpha(sub, q).diversity ≈
              norm_sub_alpha(byhand, q).diversity
    end
end

end
