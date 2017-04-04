module TestGeneralised

using Diversity
using Diversity.ρ̄

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

metacom = Metacommunity([1 2; 1 2]/6)

# Basic checks for the diversity() function
@testset "diversity()" begin
    @test diversity(metacommunityDiversity, ρ̄, metacom, 1) ≈ 1
    @test diversity(Set([subcommunityDiversity]), Set([ρ̄]), metacom, 1)[1] ≈ [1.0, 1.0]
end

end
