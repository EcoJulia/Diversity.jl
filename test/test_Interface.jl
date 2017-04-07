module TestInterface
using Base.Test

using Diversity
using Diversity.floattypes

@testset "FP" begin
    @test Float32 ∈ floattypes(Float32[1.0])
    @test Float64 ∈ floattypes(Float64[1.0])
end

end
