module TestPartition
using Test
using Diversity
using Diversity.API

three = [0.3
         0.3
         0.4]
oc_count = Onecommunity()
ab3 = [1 3 0
       2 0 4]
abnorm = ab3 / sum(ab3)
sc = Subcommunities(size(ab3, 2))
@testset "Communities" begin
    oc_2 = Onecommunity("All of it")
    @test countsubcommunities(Onecommunity("Hello")) == countsubcommunities(oc_2)
    @test getsubcommunitynames(oc_2) == ["All of it"]
    @test getsubcommunitynames(oc_count) == ["1"]
    @test getsubcommunitynames(sc) == map(x -> "$x", 1:countsubcommunities(sc))
    @test getsubcommunitynames(Subcommunities(["a", "b"])) == ["a", "b"]
end

end
