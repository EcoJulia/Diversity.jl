# SPDX-License-Identifier: BSD-2-Clause

module TestPartition
using Test
using Diversity
using Diversity.API

ab3 = [1 3 0
       2 0 4]
sc = Subcommunities(size(ab3, 2))

@testset "Communities" begin
    oc_count = Onecommunity()
    oc_2 = Onecommunity("All of it")
    @test countsubcommunities(Onecommunity("Hello")) ==
          countsubcommunities(oc_2)
    @test getsubcommunitynames(oc_2) == ["All of it"]
    @test getsubcommunitynames(oc_count) == ["1"]
    @test getsubcommunitynames(sc) == map(x -> "$x", 1:countsubcommunities(sc))
    @test getsubcommunitynames(Subcommunities(["a", "b"])) == ["a", "b"]
end

@testset "Constructor validation" begin
    # Both `Subcommunities` constructors reject an empty partition — a metacommunity divided into no
    # subcommunities has no weights to normalise by, so this must fail at construction rather than
    # producing NaNs later.
    @test_throws ErrorException Subcommunities(0)
    @test_throws ErrorException Subcommunities(-1)
    @test_throws ErrorException Subcommunities(String[])

    @test countsubcommunities(Subcommunities(1)) == 1
    @test countsubcommunities(Subcommunities(["only"])) == 1

    # `Onecommunity` is always exactly one, however it is named.
    @test countsubcommunities(Onecommunity()) == 1
    @test countsubcommunities(Onecommunity("named")) == 1
end

@testset "In a metacommunity" begin
    # A partition is only meaningful attached to abundances, and the subcommunity count has to match
    # the number of columns — the check `mcmatch` makes.
    meta = Metacommunity(ab3, UniqueTypes(2), sc)
    @test countsubcommunities(meta) == 3
    @test getsubcommunitynames(meta) == ["1", "2", "3"]
    @test_throws ErrorException Metacommunity(ab3, UniqueTypes(2),
                                              Subcommunities(2))

    # `Onecommunity` treats the whole thing as undivided, so every measure collapses to one column.
    onemeta = Metacommunity(vec(sum(ab3, dims = 2)), UniqueTypes(2),
                            Onecommunity())
    @test countsubcommunities(onemeta) == 1
    @test getweight(onemeta) ≈ [1.0]
end

end
