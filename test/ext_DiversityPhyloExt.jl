# SPDX-License-Identifier: BSD-2-Clause

module TestPhylogenetics
using Test

using Phylo
using Diversity
using Diversity.Ecology: faith_pd, generalisedfaith_pd
using EcoBase: thingkind, thingkindplural, placekind
using Tables
using DataFrames

# A phylogenetic type that is *not* this extension's `PhyloBranches` - enough of the API to build a
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

    # A types object may add columns of its own, and they have to be merged into the result before
    # it is materialised rather than inserted into a DataFrame afterwards -- a NamedTuple is
    # immutable, and a sink need not support insertion at all. So the added column has to survive
    # every route out, not just the default one.
    for f in (inddiv, subdiv, metadiv)
        columns = f(Tables.columntable, Gamma(tsmetaphylo), 0)
        @test haskey(columns, :treename)
        @test all(==("tree"), columns.treename)
        # the same number of rows as every other column, which is what makes it a valid table
        @test length(columns.treename) == length(columns.diversity)
        @test Tables.istable(typeof(f(DataFrame, Gamma(tsmetaphylo), 0)))
    end

    # And across several orders, where the columns of each are concatenated before materialising.
    many = subdiv(Gamma(tsmetaphylo), [0, 1, 2])
    @test all(==("tree"), many.treename)
    @test nrow(many) == 3 * countsubcommunities(tsmetaphylo)

    # A metacommunity whose types add nothing must not grow a column.
    plain = Metacommunity([0.4, 0.3, 0.3], ph)
    @test "treename" ∉ names(subdiv(Gamma(plain), 0))
    @test !haskey(subdiv(Tables.columntable, Gamma(plain), 0), :treename)

    # Note: The units of a phylogenetic metacommunity are *branches*, and this is how a reader is told
    # - `PhyloBranches` is opinionated about that, so it answers EcoBase's naming hooks itself.
    @test thingkind(metaphylo) == "branch"
    @test thingkindplural(metaphylo) == "branches"   # not EcoBase's default "branchs"
    @test placekind(metaphylo) == "subcommunity"
    out = sprint(show, metaphylo)
    @test occursin("with 5 branches in 1 subcommunity", out)
    @test occursin("Branch names:", out)             # singular in the heading

    # Translating to a GeneralTypes metacommunity has to carry the *scaled* Zmatrix and the
    # *branch* abundances together, or the phylogeny's numbers do not survive - which is the case
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
    # normalisation - so the numbers here are read straight off the tree, not off the framework.
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

    # It is only defined for phylogenetic types - there is no PD without a tree.
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

@testset "PhyloDistances" begin
    # Dog and Human meet 1 below the tips, Cat meets them at the root, 2 below
    species = ["Dog", "Human", "Cat"]
    nt = RootedTree(species)
    n = createnode!(nt)
    createbranch!(nt, n, species[1], 1.0)
    createbranch!(nt, n, species[2], 1.0)
    r = createnode!(nt)
    createbranch!(nt, r, n, 1.0)
    createbranch!(nt, r, species[3], 2.0)
    pd = PhyloDistances(nt)
    @test pd isa Diversity.PhyloDistances <: Diversity.AbstractPhyloTypes
    @test gettypenames(pd) == getleafnames(nt, preorder) ==
          gettypenames(PhyloBranches(nt), true)
    @test counttypes(pd) == 3
    @test getdiversityname(pd) == "Phylogenetic Distance"

    # Patristic distances 2 (Dog-Human) and 4, so the default is 1 - d / 4:
    # on an ultrametric tree, the share of the root-to-tip path in common
    order = indexin(species, gettypenames(pd))
    @test calcsimilarity(pd, 1.0)[order, order] ≈
          [1.0 0.5 0.0; 0.5 1.0 0.0; 0.0 0.0 1.0]
    named = PhyloDistances(nt, species)
    @test gettypenames(named) == species
    @test calcsimilarity(named, 1.0) ≈ [1.0 0.5 0.0; 0.5 1.0 0.0; 0.0 0.0 1.0]
    d = [0.0 2.0 4.0; 2.0 0.0 4.0; 4.0 4.0 0.0]
    @test calcsimilarity(PhyloDistances(nt, species; transform = :exponential),
                         1.0) ≈ exp.(-d ./ 4)
    @test calcsimilarity(PhyloDistances(nt, species; k = 0.1,
                                        normalise = false),
                         1.0) ≈ 1 .- 0.1 .* d
    @test calcsimilarity(PhyloDistances(nt, species; k = 0.5,
                                        normalise = false),
                         1.0) ≈ max.(1 .- 0.5 .* d, 0)

    # Species are the types, so diversity is that of GeneralTypes with Z
    Z = calcsimilarity(named, 1.0)
    abund = [0.25 0.05; 0.15 0.05; 0.1 0.4]
    mc = Metacommunity(abund, named)
    byhand = Metacommunity(abund, GeneralTypes(Z, species))
    @test counttypes(mc) == 3
    @test getordinariness!(mc) ≈ Z * abund
    for q in [0, 1, 2, Inf]
        @test meta_gamma(mc, q).diversity ≈ meta_gamma(byhand, q).diversity
        @test norm_sub_rho(mc, q).diversity ≈ norm_sub_rho(byhand, q).diversity
    end
    # and a subset of species keeps the similarity as it was
    sub = view(mc, species = [1, 3])
    @test calcsimilarity(gettypes(sub), 1.0) ≈ Z[[1, 3], [1, 3]]

    # Any rooted tree: 1 - d / max d from the tree's own distances
    for tree in (rand(Nonultrametric(20)), rand(Ultrametric(20)))
        pdt = PhyloDistances(tree)
        names = gettypenames(pdt)
        dt = [distance(tree, a, b) for a in names, b in names]
        @test calcsimilarity(pdt, 1.0) ≈ 1 .- dt ./ maximum(dt)
    end

    # From a TreeSet of one tree, with the tree named in the output
    ts = PhyloDistances(TreeSet(Dict("tree" => nt)), species)
    @test calcsimilarity(ts, 1.0) == Z
    @test meta_gamma(Metacommunity(abund, ts), 1).treename == ["tree"]

    # Every leaf exactly once, as for PhyloBranches
    @test_throws ArgumentError PhyloDistances(nt, ["Dog", "Cat"])
    @test_throws ArgumentError PhyloDistances(nt, ["Dog", "Dog", "Cat"])
    @test_throws ArgumentError PhyloDistances(nt, [species; "Mouse"])
    @test_throws ArgumentError PhyloDistances(TreeSet(Dict("tree" => nt)),
                                              ["Dog", "Cat"])
    @test_throws ArgumentError PhyloDistances(nt; transform = :cubic)
    nolengths = RootedTree(species)
    m = createnode!(nolengths)
    createbranch!(nolengths, m, "Dog")
    createbranch!(nolengths, m, "Human")
    top = createnode!(nolengths)
    createbranch!(nolengths, top, m)
    createbranch!(nolengths, top, "Cat")
    @test_throws "needs a length" PhyloDistances(nolengths)
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

    # The subset measures with its *own* scale, not the parent's. This is the one place a view has
    # to compute something rather than carry it: the scale is the abundance-weighted mean
    # root-to-tip distance, so dropping subcommunities changes it, and it cannot be recovered from
    # branch abundances alone. Getting it wrong is not subtle -- inheriting the generic default of
    # one doubled every diversity here.
    @test meta_gamma(sites, 1).diversity ≈
          meta_gamma(Metacommunity(getabundance(mc, true)[:, 1:1] ./
                                   sum(getabundance(mc, true)[:, 1:1]), ph),
                     1).diversity

    # And on a non-ultrametric tree the parent's scale and the subset's genuinely differ, so this
    # would fail if the parent's were carried across instead of the subset's being worked out.
    nu = RootedTree(species)
    ni = createnode!(nu)
    createbranch!(nu, ni, species[1], 1.0)
    createbranch!(nu, ni, species[2], 3.0)
    nr = createnode!(nu)
    createbranch!(nu, nr, ni, 1.0)
    createbranch!(nu, nr, species[3], 2.0)
    nuleaf = [0.5 0.02; 0.02 0.02; 0.02 0.42]
    numc = Metacommunity(nuleaf, PhyloBranches(nu))
    nuview = view(numc, sites = 1:1)
    @test Diversity.API._getscale(nuview) ≉ Diversity.API._getscale(numc)
    @test meta_gamma(nuview, 1).diversity ≈
          meta_gamma(Metacommunity(nuleaf[:, 1:1] ./ sum(nuleaf[:, 1:1]),
                                   PhyloBranches(nu)), 1).diversity

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

@testset "leaf order" begin
    # Non-ultrametric, so every branch length matters, and its preorder (a, b, d, e, c) is not even
    # the order the Newick string writes its leaves in.
    tree = parsenewick("((a:1,b:2):1,(c:1,(d:0.5,e:3):1):2);")
    ph = PhyloBranches(tree)
    order = ["e", "a", "d", "b", "c"]
    perm = [findfirst(==(leaf), gettypenames(ph, true)) for leaf in order]
    reordered = PhyloBranches(tree, order)

    # The species take the caller's order, and the branches follow them, grouped by species.
    @test gettypenames(reordered, true) == order
    @test sort(gettypenames(reordered)) == sort(gettypenames(ph))
    owners = [first(split(name, " : ")) for name in gettypenames(reordered)]
    @test issorted(owners, by = leaf -> findfirst(==(leaf), order))

    # One community with its rows in either order measures the same: the rows here are in the
    # tree's own order, and each subcommunity lacks a different species.
    abund = [0.3 0.0 0.1;
             0.1 0.2 0.0;
             0.2 0.1 0.0;
             0.1 0.0 0.3;
             0.0 0.3 0.2]
    abund ./= sum(abund)
    mc = Metacommunity(abund, ph)
    mcre = Metacommunity(abund[perm, :], reordered)
    for q in [0, 1, 2, Inf]
        @test meta_gamma(mcre, q).diversity ≈ meta_gamma(mc, q).diversity
        @test norm_sub_alpha(mcre, q).diversity ≈
              norm_sub_alpha(mc, q).diversity
        @test norm_sub_beta(mcre, q).diversity ≈ norm_sub_beta(mc, q).diversity
    end
    @test faith_pd(mcre).diversity ≈ faith_pd(mc).diversity
    # The same rows read in the tree's order are a different community.
    @test !(faith_pd(Metacommunity(abund[perm, :], ph)).diversity ≈
            faith_pd(mc).diversity)

    # Naming the tree's own order is the default.
    same = PhyloBranches(tree, gettypenames(ph, true))
    @test gettypenames(same) == gettypenames(ph)
    @test calcsimilarity(same, 1.0) == calcsimilarity(ph, 1.0)

    tsre = PhyloBranches(TreeSet(Dict("tree" => tree)), order)
    @test gettypenames(tsre, true) == order
    @test gettypenames(tsre) == gettypenames(reordered)
    @test calcsimilarity(tsre, 1.0) == calcsimilarity(reordered, 1.0)
    tsmc = Metacommunity(abund[perm, :], tsre)
    @test meta_gamma(tsmc, 1).diversity ≈ meta_gamma(mc, 1).diversity
    @test meta_gamma(tsmc, 1).treename == ["tree"]

    # Every leaf exactly once, and nothing else.
    @test_throws ArgumentError PhyloBranches(tree, ["e", "a", "d", "b"])
    @test_throws ArgumentError PhyloBranches(tree, [order; "f"])
    @test_throws ArgumentError PhyloBranches(tree, [order; "a"])
    @test_throws ArgumentError PhyloBranches(tree, ["e", "a", "d", "b", "b"])
    internal = getnodename(tree, getparent(tree, "a"))
    @test_throws ArgumentError PhyloBranches(tree, [order; internal])
    @test_throws ArgumentError PhyloBranches(TreeSet(Dict("tree" => tree)),
                                             ["e", "a", "d", "b"])
    @test gettypenames(PhyloBranches(tree, split("e a d b c")), true) == order

    # PhyloDistances takes its species order the same way, and by default the same order as
    # PhyloBranches, so one abundance matrix means the same community under either type.
    pd = PhyloDistances(tree)
    @test gettypenames(pd) == gettypenames(ph, true)
    pdre = PhyloDistances(tree, order)
    @test gettypenames(pdre) == order
    @test calcsimilarity(pdre, 1.0) ≈ calcsimilarity(pd, 1.0)[perm, perm]
    mcd = Metacommunity(abund, pd)
    mcdre = Metacommunity(abund[perm, :], pdre)
    for q in [0, 1, 2, Inf]
        @test meta_gamma(mcdre, q).diversity ≈ meta_gamma(mcd, q).diversity
        @test norm_sub_beta(mcdre, q).diversity ≈
              norm_sub_beta(mcd, q).diversity
    end
    tsd = PhyloDistances(TreeSet(Dict("tree" => tree)), order)
    @test gettypenames(tsd) == order
    @test calcsimilarity(tsd, 1.0) == calcsimilarity(pdre, 1.0)
    @test gettypenames(PhyloDistances(tree, split("e a d b c"))) == order
end

end
