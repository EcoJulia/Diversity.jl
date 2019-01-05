module TestEffectiveNumbers
using Compat.Test
using Compat.LinearAlgebra
using Compat

using Diversity
using Diversity: powermean

numbers = [1.0, 2, 4, 8, 16];
numspecies = 100;
fragments = rand(numspecies);
weights = rand(numspecies);
weights /= sum(weights);
Z1 = ones(typeof(weights[1]), (length(weights), length(weights)));
numcommunities = 8;
manyweights = rand(numspecies, numcommunities);
manyweights *= Diagonal(reshape(mapslices(v -> 1. / sum(v), manyweights, dims=1),
                                numcommunities));

# Simple power means - we no longer export these, but we should check
# them anyway as everything relies on them

@testset "powermean" begin
    # Check that an exception is thrown when 'values' and weights are different lengths
    @test_throws DimensionMismatch powermean(numbers, 0, weights)

    # Some simple values
    @test powermean([1.0], 0.0, [1.0]) ≈ 1.0
    @test powermean(numbers, 0.0) ≈ 4.0
    @test powermean(numbers, [-Inf]) ≈ [1]
    @test powermean(numbers, [1.0, -1.0]) ≈ [31.0/5.0, 80.0/31.0]
    @test powermean(numbers, Inf, [1.0, 1.0, 1.0, 1.0, 0.0]) ≈ 8
    @test isnan(powermean(numbers, 0.0, 0.0 * numbers))

    # Power mean with some random numbers
    @test powermean(fragments, 0) ≈ prod(fragments .^ (1. / numspecies))
    @test powermean(fragments, 1) ≈ Compat.Statistics.mean(fragments)
    @test powermean(fragments, Inf) ≈ maximum(fragments)
    @test powermean(fragments, 0, weights) ≈ prod(fragments .^ weights)
    @test powermean(fragments, 1, weights) ≈ sum(fragments .* weights)
    @test powermean(manyweights, -1, manyweights) .^ -1 ≈ numspecies * ones(size(manyweights, 2))
end

@testset "qD" begin
    # Basic qD diversity calculation
    @test qD(weights, 0) ≈ mapreduce((x) -> x ≈ 0 ? 0 : 1, +, weights)
    @test qD(Metacommunity(weights), 0) == qD(weights, 0)
    @test qD(Metacommunity(weights, UniqueTypes(numspecies)), 0) ==
        qD(weights, 0)
    @test_throws ErrorException qD(Metacommunity(weights,
                                   GeneralTypes(rand(numspecies, numspecies))),
                                   0)
    @test qD(weights, 1) ≈ prod(weights .^ -weights)
    @test qD(weights, 2) ≈ 1.0 / sum(weights .^ 2)
    @test qD(weights, Inf) ≈ 1.0 / maximum(weights)

    @test qD(weights, [1, 2]) ≈ [qD(weights, 1), qD(weights, 2)]

    @test typeof(qD(manyweights[:,1], 0)) <: AbstractFloat
    @test typeof(qD(manyweights[:,1], [0])) <: Vector

    for i in 1:size(manyweights, 2)
        @test qD(manyweights[:,i], [0]) ≈ numspecies * ones((1, size(manyweights[:,i], 2)))
    end

    # Diversities are not normalised, so generate an error
    if VERSION < v"0.7.0-"
        @test_warn "Abundances not normalised to 1, correcting..." qD([0.1, 0.1], 1)
    end
end

@testset "qDZ" begin
    # General Leinster-Cobbold diversity calculation
    @test qDZ(weights, [1, 2]) ≈ qD(weights, [1, 2])
    @test qDZ(weights, [0, 1, 2, 3, Inf], Z1) ≈ [1, 1, 1, 1, 1]


    for i in 1:size(manyweights, 2)
        @test qDZ(manyweights[:,i], [0, 1, 2, Inf],
                  ones((size(manyweights[:,i], 1),
                        size(manyweights[:,i], 1)))) ≈ ones((4, size(manyweights[:,i], 2)))
    end

    if VERSION < v"0.7.0-"
        # Diversities are not normalised, so generate an warning
        @test_warn "Abundances not normalised to 1, correcting..." qDZ([0.1, 0.1], 1)
    end
end

end
