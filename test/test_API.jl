module TestAPI
using Compat.Test

using Diversity.API

@testset "FP" begin
    @test Float32 ∈ floattypes(Float32[1.0])
    @test Float64 ∈ floattypes(Float64[1.0])
end

end
