using Diversity
using Base.Test

# Simple power means
numbers = [1, 2, 4, 8, 16];
@test_approx_eq powermean(numbers, 0) 4
@test_approx_eq powermean(numbers, 1) 6.2
@test_approx_eq powermean(numbers, -Inf) 1
@test_approx_eq powermean(numbers, Inf, [1, 1, 1, 1, 0]) 8

# Power mean with some random numbers
len = 100;
fragments = rand(len);
weights = rand(len);
weights = weights / sum(weights);
@test_approx_eq powermean(fragments, 0) prod(fragments .^ (1/len))
@test_approx_eq powermean(fragments, 1) mean(fragments)
@test_approx_eq powermean(fragments, Inf) maximum(fragments)
@test_approx_eq powermean(fragments, 0, weights) prod(fragments .^ weights)
@test_approx_eq powermean(fragments, 1, weights) sum(fragments .* weights)

# Basic qD diversity calculation
@test_approx_eq qD(weights, 0) mapreduce((x) -> isapprox(x, 0.) ? 0. : 1.,
                                         +, weights)
@test_approx_eq qD(weights, 1) prod(weights .^ -weights)
@test_approx_eq qD(weights, 2) 1. / sum(weights .^ 2)
@test_approx_eq qD(weights, Inf) 1. / maximum(weights)

@test_approx_eq qD(weights, [1, 2]) [qD(weights, 1), qD(weights, 2)]
