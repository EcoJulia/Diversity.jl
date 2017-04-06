module TestGeneralised

using Diversity
using Diversity.ρ̄
using Base.Test

metacom = Metacommunity([1 2; 1 2]/6)

# Basic checks for the diversity() function
@testset "diversity()" begin
    @test diversity([metacommunityDiversity], [ρ̄], metacom, 1)[:diversity] ≈ [1.0]
    @test diversity(Set([subcommunityDiversity]), Set([ρ̄]), metacom, 1)[:diversity] ≈ [1.0, 1.0]
end

end
