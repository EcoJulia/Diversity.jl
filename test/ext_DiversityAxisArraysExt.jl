module TestAxisArrays
using Test

using AxisArrays
using Diversity
using LinearAlgebra

@testset "AxisArrays" begin
    species = ["Dog", "Human", "Cat"]
    A = AxisArray(I(3) .* 1.0, Axis{:row}(species), Axis{:col}(species))
    a = GeneralTypes(A)
    @test gettypenames(a) == species
    B = AxisArray(I(3) .* -1.0, Axis{:row}(species), Axis{:col}(species))
    @test_throws "Similarities must be ≥ 0" GeneralTypes(B)
    C = AxisArray(I(3) .* 1.0, Axis{:row}(species), Axis{:col}(sort(species)))
    @test_throws DimensionMismatch GeneralTypes(C)
end

end
