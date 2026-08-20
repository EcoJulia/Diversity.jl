# SPDX-License-Identifier: BSD-2-Clause

module TestEffectiveNumbers
using Test
using LinearAlgebra
using Statistics

using Diversity
using Diversity: powermean, _powermeancolumns, THREADINGTHRESHOLD
using Random

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

@testset "Threaded power means" begin
    # `powermean` over a matrix reduces each column independently, so above a size threshold the
    # columns are spread over threads. The two paths must agree exactly -- this is a reduction with
    # no randomness, so thread count cannot change the answer, and what the threaded path can get
    # wrong is the order it writes results back in.
    #
    # Warning: the test workers run with JULIA_NUM_THREADS=1, so `powermean` itself would never
    # choose the threaded path here. Both are called directly for that reason.
    Random.seed!(1)
    for (nt, np) in ((1, 1), (3, 1), (1, 5), (10, 7), (25, 40))
        values = rand(nt, np) .+ 0.5
        weights = rand(nt, np)
        for order in (0, 0.5, 1, 2, -1, Inf, -Inf)
            serial = _powermeancolumns(values, order, weights, false)
            threaded = _powermeancolumns(values, order, weights, true)
            @test serial == threaded
            @test length(threaded) == np
        end
        # A vector of orders gives a vector per column, so the threaded path has to get the
        # element type right as well as the order.
        orders = [0, 1, 2]
        @test _powermeancolumns(values, orders, weights, false) ==
              _powermeancolumns(values, orders, weights, true)
        @test eltype(_powermeancolumns(values, orders, weights, true)) <:
              AbstractVector
    end

    # A column of zero weights is NaN by definition, and NaN != NaN, so the paths are compared
    # elementwise here rather than by equality.
    values = rand(4, 3) .+ 0.5
    weights = rand(4, 3)
    weights[:, 2] .= 0.0
    serial = _powermeancolumns(values, 1, weights, false)
    threaded = _powermeancolumns(values, 1, weights, true)
    @test isnan(threaded[2])
    @test all(isequal(s, t) for (s, t) in zip(serial, threaded))

    # And through the public function, on a matrix big enough to cross the threshold, so that
    # whichever path it picks it still agrees with the serial one.
    big = rand(20, cld(THREADINGTHRESHOLD, 20) + 10) .+ 0.5
    bigw = rand(size(big)...)
    @test length(big) > THREADINGTHRESHOLD
    @test powermean(big, 1, bigw) ==
          _powermeancolumns(big, 1, bigw, false)
end

end
