# SPDX-License-Identifier: BSD-2-Clause

module TestInterface
using Test

using Diversity
using Diversity: createsummaryline
using EcoBase
using EcoBase: thingkind, thingkindplural, placekind, placekindplural
using LinearAlgebra

numspecies = 10
numcommunities = 8
manyweights = rand(numspecies, numcommunities)
manyweights /= sum(manyweights)

@testset "createsummaryline" begin
    # One name is itself; a short list is comma-separated; a long one elides the middle. All three
    # branches, because only the middle one was ever exercised.
    @test createsummaryline(["only"]) == "only"
    @test createsummaryline(["a", "b"]) == "a, b"
    @test createsummaryline(["a", "b", "c", "d", "e"]) == "a, b, c, d, e"

    long = createsummaryline(["a", "b", "c", "d", "e", "f", "g"])
    @test occursin("...", long)
    @test startswith(long, "a, b, c")
    @test endswith(long, "f, g")

    # Type names are not necessarily strings — `GeneralTypes(zmatrix)` numbers its types from the
    # matrix axes — and showing such a metacommunity threw a `MethodError` until this worked.
    @test createsummaryline([1, 2, 3]) == "1, 2, 3"
end

@testset "Text output" begin
    species = map(n -> "Species $n", 1:numspecies)
    communities = map(n -> "SC $n", 1:numcommunities)
    ut = UniqueTypes(species)
    sc = Subcommunities(communities)
    mc = Metacommunity(manyweights, ut, sc)

    io = IOBuffer()
    show(io, mc)
    out = String(take!(io))
    @test occursin("measuring", out)
    @test occursin("Unique", out)
    @test occursin("Species 1", out)
    @test occursin("SC 1", out)

    # The unnamed case, which is what a bare `Metacommunity(abundances, Z)` gives you.
    io = IOBuffer()
    show(io, Metacommunity(manyweights, Matrix(1.0I, numspecies, numspecies)))
    @test occursin("Arbitrary Z", String(take!(io)))
end

@testset "Accessors" begin
    species = map(n -> "Species $n", 1:numspecies)
    communities = map(n -> "SC $n", 1:numcommunities)
    mc = Metacommunity(manyweights, UniqueTypes(species),
                       Subcommunities(communities))

    # Asserted directly rather than left to incidental coverage from other files — these are the
    # package's public reading surface, and a file that happened to exercise them could move.
    @test counttypes(mc) == numspecies
    @test countsubcommunities(mc) == numcommunities
    @test gettypenames(mc) == species
    @test getsubcommunitynames(mc) == communities
    @test counttypes(gettypes(mc)) == numspecies
    @test countsubcommunities(getpartition(mc)) == numcommunities
    @test getdiversityname(mc) == "Unique"
    @test !hassimilarity(mc)

    # Abundances are relative to the whole metacommunity; weights are the subcommunity column sums;
    # the metacommunity abundance is the row sums. All three must agree with each other.
    @test sum(getabundance(mc)) ≈ 1.0
    @test getweight(mc) ≈ vec(sum(manyweights, dims = 1))
    @test getmetaabundance(mc) ≈ vec(sum(manyweights, dims = 2))
    @test sum(getweight(mc)) ≈ 1.0

    # With no similarity, ordinariness is the abundance itself.
    @test getordinariness!(mc) ≈ getabundance(mc)
    @test getmetaordinariness!(mc) ≈ getmetaabundance(mc)

    # `raw` selects the abundances as supplied rather than the processed ones. They coincide here,
    # because nothing rescales them, but the argument must still be honoured.
    @test getabundance(mc, true) ≈ manyweights
    @test gettypenames(mc, true) == gettypenames(mc, false)
    @test counttypes(mc, true) == counttypes(mc, false)

    # No added output columns unless a type asks for them (the Phylo extension does).
    @test isempty(addedoutputcols(mc))
    @test isnothing(getaddedoutput(mc))
end

@testset "Naming the units" begin
    # EcoBase asks the assemblage what its units are called and defaults to "thing"/"place"; we
    # answer on the types and the partition instead, so that output says what it means. This is
    # also how a reader is told that a phylogeny's units are branches rather than species.
    mc = Metacommunity(manyweights)
    @test thingkind(mc) == "species"
    @test placekind(mc) == "subcommunity"
    @test placekindplural(mc) == "subcommunities"
    @test thingkindplural(mc) == "species"

    # It is the *types* that are asked, not the metacommunity, so a new types object can answer.
    @test thingkind(gettypes(mc)) == thingkind(mc)
    @test placekind(getpartition(mc)) == placekind(mc)

    out = sprint(show, mc)
    @test occursin("with $numspecies species in $numcommunities subcommunities",
                   out)
    @test occursin("Species names:", out)
    @test occursin("Subcommunity names:", out)
    @test occursin("measuring Unique diversity", out)

    # And the singular really is used when there is one of something.
    single = sprint(show, Metacommunity([1.0]))
    @test occursin("with 1 species in 1 subcommunity measuring", single)
end

end
