# SPDX-License-Identifier: BSD-2-Clause

module TestEffectiveNumbers
using Test
using LinearAlgebra
using Statistics

using Diversity
using Diversity: powermean

numbers = [1.0, 2, 4, 8, 16]
numspecies = 100
fragments = rand(numspecies)
weights = rand(numspecies)
weights /= sum(weights)
Z1 = ones(typeof(weights[1]), (length(weights), length(weights)))
numcommunities = 8
manyweights = rand(numspecies, numcommunities)
manyweights *= Diagonal(reshape(mapslices(v -> 1.0 / sum(v), manyweights,
                                          dims = 1),
                                numcommunities))

# Simple power means - we no longer export these, but we should check
# them anyway as everything relies on them

@testset "powermean" begin
    # Check that an exception is thrown when 'values' and weights are different lengths
    @test_throws DimensionMismatch powermean(numbers, 0, weights)

    # Some simple values
    @test powermean([1.0], 0.0, [1.0]) ≈ 1.0
    @test powermean(numbers, 0.0) ≈ 4.0
    @test powermean(numbers, [-Inf]) ≈ [1]
    @test powermean(numbers, [1.0, -1.0]) ≈ [31.0 / 5.0, 80.0 / 31.0]
    @test powermean(numbers, Inf, [1.0, 1.0, 1.0, 1.0, 0.0]) ≈ 8
    @test isnan(powermean(numbers, 0.0, 0.0 * numbers))

    # Power mean with some random numbers
    @test powermean(fragments, 0) ≈ prod(fragments .^ (1.0 / numspecies))
    @test powermean(fragments, 1) ≈ mean(fragments)
    @test powermean(fragments, Inf) ≈ maximum(fragments)
    @test powermean(fragments, 0, weights) ≈ prod(fragments .^ weights)
    @test powermean(fragments, 1, weights) ≈ sum(fragments .* weights)
    @test powermean(manyweights, -1, manyweights) .^ -1 ≈
          numspecies * ones(size(manyweights, 2))
end

@testset "qD" begin
    # Basic qD diversity calculation
    @test qD(weights, 0) ≈ mapreduce((x) -> x ≈ 0 ? 0 : 1, +, weights)
    @test qD(Metacommunity(weights), 0) == qD(weights, 0)
    @test qD(Metacommunity(weights, UniqueTypes(numspecies)), 0) ==
          qD(weights, 0)
    @test_throws ErrorException qD(Metacommunity(weights,
                                                 GeneralTypes(rand(numspecies,
                                                                   numspecies))),
                                   0)
    @test qD(weights, 1) ≈ prod(weights .^ -weights)
    @test qD(weights, 2) ≈ 1.0 / sum(weights .^ 2)
    @test qD(weights, Inf) ≈ 1.0 / maximum(weights)

    @test qD(weights, [1, 2]) ≈ [qD(weights, 1), qD(weights, 2)]

    @test typeof(qD(manyweights[:, 1], 0)) <: AbstractFloat
    @test typeof(qD(manyweights[:, 1], [0])) <: Vector

    for i in axes(manyweights, 2)
        @test qD(manyweights[:, i], [0]) ≈
              numspecies * ones((1, size(manyweights[:, i], 2)))
    end
end

@testset "qDZ" begin
    # General Leinster-Cobbold diversity calculation
    @test qDZ(weights, [1, 2]) ≈ qD(weights, [1, 2])
    @test qDZ(weights, [0, 1, 2, 3, Inf], Z1) ≈ [1, 1, 1, 1, 1]

    for i in axes(manyweights, 2)
        @test qDZ(manyweights[:, i], [0, 1, 2, Inf],
                  ones((size(manyweights[:, i], 1),
                        size(manyweights[:, i], 1)))) ≈
              ones((4, size(manyweights[:, i], 2)))
    end
end

@testset "Empty subcommunities are recognised from their weights" begin
    # A subcommunity of zero weight holds no individuals, so its power mean is NaN. The vector
    # method discovers that by scanning the whole column; handing the matrix method the
    # per-subcommunity weights lets it read the same fact off a number it already had. The two must
    # agree exactly -- including the NaNs, hence `isequal` rather than `==`.
    values = rand(6, 5) .+ 0.5
    weights = rand(6, 5)
    weights[:, 2] .= 0.0
    weights[:, 5] .= 0.0
    colweights = vec(sum(weights, dims = 1))
    @test iszero(colweights[2]) && iszero(colweights[5])

    for order in (0, 0.5, 1, 2, -1, Inf, -Inf)
        scanned = powermean(values, order, weights)
        told = powermean(values, order, weights, colweights)
        @test isequal(scanned, told)
        @test isnan(told[2]) && isnan(told[5])
        @test !isnan(told[1])
    end

    # A vector of orders gives a vector per subcommunity, so an empty one has to come back the
    # right shape as well as the right value.
    orders = [0, 1, 2]
    @test isequal(powermean(values, orders, weights),
                  powermean(values, orders, weights, colweights))
    told = powermean(values, orders, weights, colweights)
    @test length(told[2]) == length(orders)
    @test all(isnan, told[2])

    # One weight per subcommunity, or say so here rather than failing on an index somewhere later.
    @test_throws DimensionMismatch powermean(values, 1, weights,
                                             colweights[1:3])
end

end
