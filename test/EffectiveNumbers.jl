module TestEffectiveNumbers
using Diversity
using Compat
using Base.Test

numbers = [1.0, 2, 4, 8, 16];
numspecies = 100;
fragments = rand(numspecies);
weights = rand(numspecies);
weights /= sum(weights);
Z1 = ones(typeof(weights[1]), (length(weights), length(weights)));
numcommunities = 8;
manyweights = rand(numspecies, numcommunities);
manyweights *= diagm(reshape(mapslices(v -> 1. / sum(v), manyweights, 1),
                             (numcommunities)));

# Simple power means - we no longer export these, but we should check
# them anyway as everything relies on them

# Check that an exception is thrown when 'values' and weights are different lengths
@test_throws DimensionMismatch Diversity.powermean(numbers, 0, weights)

# Some simple values
@test_approx_eq Diversity.powermean(numbers, 0.0) 4.0
@test_approx_eq Diversity.powermean(numbers, [-Inf]) [1]
@test_approx_eq Diversity.powermean(numbers, [1.0, -1.0]) [31/5, 80/31]
@test_approx_eq Diversity.powermean(numbers, Inf, [1.0, 1.0, 1.0, 1.0, 0.0]) 8
@test isnan(Diversity.powermean(numbers, 0.0, 0.0 * numbers))

# Power mean with some random numbers
@test_approx_eq Diversity.powermean(fragments, 0) prod(fragments .^ (1. / numspecies))
@test_approx_eq Diversity.powermean(fragments, 1) mean(fragments)
@test_approx_eq Diversity.powermean(fragments, Inf) maximum(fragments)
@test_approx_eq Diversity.powermean(fragments, 0, weights) prod(fragments .^ weights)
@test_approx_eq Diversity.powermean(fragments, 1, weights) sum(fragments .* weights)

# Basic qD diversity calculation
@test_approx_eq qD(weights, 0) mapreduce((x) -> isapprox(x, 0) ? 0 : 1, +, weights)
@test_approx_eq qD(weights, 1) prod(weights .^ -weights)
@test_approx_eq qD(weights, 2) 1. / sum(weights .^ 2)
@test_approx_eq qD(weights, Inf) 1. / maximum(weights)

@test_approx_eq qD(weights, [1, 2]) [qD(weights, 1), qD(weights, 2)]

# General Leinster-Cobbold diversity calculation
@test_approx_eq qDZ(weights, [1, 2]) qD(weights, [1, 2])
@test_approx_eq qDZ(weights, [0, 1, 2, 3, Inf], Z1) [1, 1, 1, 1, 1]

@test typeof(qD(manyweights[:,1], 0)) <: AbstractFloat
@test typeof(qD(manyweights[:,1], [0])) <: Vector

for i in 1:size(manyweights, 2)
    @test_approx_eq qD(manyweights[:,i], [0]) numspecies * ones((1, size(manyweights[:,i], 2)))
    @test_approx_eq qDZ(manyweights[:,i], [0, 1, 2, Inf],
                        ones((size(manyweights[:,i], 1),
                              size(manyweights[:,i], 1)))) ones((4, size(manyweights[:,i], 2)))
end

# Generate warnings, but normalise and calculate diversities
warn("We now generate two warnings for code coverage completeness...")
@test_approx_eq qD([0.1, 0.1], 1) 2.
@test_approx_eq qDZ([0.1, 0.1], 1) 2.

end

