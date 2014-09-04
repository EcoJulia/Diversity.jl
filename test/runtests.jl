using Diversity
using Base.Test

# Simple power means - we no longer export these, but we should check
# them anyway as everything relies on them
using Diversity.powermean
numbers = [1, 2, 4, 8, 16];
@test_approx_eq powermean(numbers, 0) 4
@test_approx_eq powermean(numbers, 1) 6.2
@test_approx_eq powermean(numbers, -Inf) 1
@test_approx_eq powermean(numbers, Inf, [1, 1, 1, 1, 0]) 8

# Power mean with some random numbers
numspecies = 100;
fragments = rand(numspecies);
weights = rand(numspecies);
weights /= sum(weights);
@test_throws ErrorException powermean(numbers, 0, weights)
@test_approx_eq powermean(fragments, 0) prod(fragments .^ (1. / numspecies))
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

# General Leinster-Cobbold diversity calculation
Z1 = ones(typeof(weights[1]), (length(weights), length(weights)));
@test_approx_eq qDZ(weights, [1, 2]) qD(weights, [1, 2])
@test_approx_eq qDZ(weights, [0, 1, 2, 3, Inf], Z1) [1, 1, 1, 1, 1]

numcommunities = 8;
manyweights = rand(numspecies, numcommunities);
manyweights *= diagm(reshape(mapslices(v -> 1. / sum(v), manyweights, 1),
                             (numcommunities)));

@test_approx_eq qD(manyweights, 0) numspecies * ones((1, size(manyweights)[2]))
@test_approx_eq qD(manyweights, [0]) numspecies * ones((1, size(manyweights)[2]))
@test_approx_eq qDZ(manyweights, [0, 1, 2, Inf],
                    ones((size(manyweights)[1],
                          size(manyweights)[1]))) ones((4, size(manyweights)[2]))

# Generate warnings, but normalise and calculate diversities
warn("We now generate two warnings for code coverage completeness...")
@test_approx_eq qD([0.1, 0.1], 1) 2.
@test_approx_eq qDZ([0.1, 0.1], 1) 2.

# Sub-community alpha diversities
communities = rand(numspecies, numcommunities);
communities /= sum(communities);
@test_approx_eq ᾱ(communities, 0) numspecies * ones((1, size(communities)[2]))
@test_approx_eq communityalphabar(communities,
                                  [0]) numspecies * ones((1, size(communities)[2]))
@test_approx_eq ᾱ(communities,
                  [0, 1, 2, Inf], Z1) ones((4, size(communities)[2]))

@test_approx_eq α(communities, 0) numspecies * mapslices(v -> 1. / sum(v),
                                                  communities, 1)

even = ones((numspecies, numcommunities)) / (numspecies * numcommunities);
qs = [0, 1, 2, 3, 4, 5, 6, Inf];
@test_approx_eq Ā(even, qs) numspecies * ones((1, length(qs)))
@test_approx_eq A(even, qs) numspecies * numcommunities * ones((1, length(qs)))

probs = reshape(mapslices(sum, communities, 2), (size(communities)[1]));
@test_approx_eq G(communities, qs) Ḡ(communities, qs)
@test_approx_eq G(communities, qs) qD(probs, qs)
@test_approx_eq G(communities, qs, Z1) qDZ(probs, qs, Z1)

Z = rand(numspecies, numspecies);
@test_approx_eq G(communities, qs, Z) qDZ(probs, qs, Z)

colweights = rand(numcommunities);
colweights /= sum(colweights);
allthesame = probs * colweights';
@test_approx_eq B̄(allthesame, qs, Z) ones((1, length(qs)))
@test_approx_eq B(allthesame, qs) powermean(colweights, 1 - qs, colweights)

# Need to check the diversity() function
@test_approx_eq diversity(B̄, allthesame, qs, Z, false, false, true) colweights
@test_approx_eq diversity(β̄, communities, qs, Z,
                          false, true, false) β̄(communities, qs, Z)
@test_approx_eq diversity(β̄, communities, qs, Z,
                          true, false, false) B̄(communities, qs, Z)
@test_approx_eq diversity(ᾱ, communities, qs, Z,
                          true, true, false)[1] Ā(communities, qs, Z)
@test_approx_eq diversity(α, allthesame, qs, Z,
                          true, true, true)[2] α(allthesame, qs, Z)
ed, cd, w = diversity(γ, communities, qs, Z, true, true, true)
@test_approx_eq diversity(γ̄, communities, qs, Z, true, false, true)[1] ed 
@test_approx_eq diversity(γ̄, communities, qs, Z, false, true, true)[1] cd
@test_approx_eq diversity(α, allthesame, qs, Z, true, true, true)[3] colweights

# Now some even communities, should see that raw and normalised
# diversities are the same
smoothed = communities ./ mapslices(sum, communities, 1);
smoothed /= numcommunities;
# Just for completeness, check one for q=-Inf - we currently have no use for this, but it is coded.
@test_approx_eq contributions(α, smoothed, [-Inf, 0:5, Inf], true) contributions(ᾱ, smoothed, [-Inf, 0:5, Inf], true)
@test_approx_eq contributions(β, smoothed, [0:5, Inf], true) contributions(β̄, smoothed, [0:5, Inf], true)
@test_approx_eq contributions(γ, smoothed, [0:5, Inf], true) contributions(γ̄, smoothed, [0:5, Inf], true)
@test_approx_eq contributions(α, smoothed, [0:5, Inf], false) contributions(ᾱ, smoothed, [0:5, Inf], false)
@test_approx_eq contributions(β, smoothed, [0:5, Inf], false) contributions(β̄, smoothed, [0:5, Inf], false)
@test_approx_eq contributions(γ, smoothed, [0:5, Inf], false) contributions(γ̄, smoothed, [0:5, Inf], false)

@test_approx_eq contributions(α, smoothed, [0:5, Inf], true) contributions(α, smoothed, [0:5, Inf], false) * numcommunities
@test_approx_eq contributions(β, smoothed, [0:5, Inf], true) contributions(β, smoothed, [0:5, Inf], false) * numcommunities
@test_approx_eq contributions(γ, smoothed, [0:5, Inf], true) contributions(γ, smoothed, [0:5, Inf], false) * numcommunities
@test_approx_eq contributions(ᾱ, smoothed, [0:5, Inf], true) contributions(ᾱ, smoothed, [0:5, Inf], false) * numcommunities
@test_approx_eq contributions(β̄, smoothed, [0:5, Inf], true) contributions(β̄, smoothed, [0:5, Inf], false) * numcommunities
@test_approx_eq contributions(γ̄, smoothed, [0:5, Inf], true) contributions(γ̄, smoothed, [0:5, Inf], false) * numcommunities

# Looking at relations to historical measures, updated with similarity
# and partitioning
using Diversity.Compatibility

@test_approx_eq richness(communities) ᾱ(communities, 0)
@test_approx_eq generalisedrichness(Ḡ, communities, Z1) 1

@test_approx_eq shannon(communities) log(ᾱ(communities, 1))
@test_approx_eq generalisedshannon(Ḡ, communities, Z1) 0

@test_approx_eq simpson(communities) ᾱ(communities, 2) .^ -1
@test_approx_eq generalisedsimpson(Ḡ, communities, Z1) 1

# Checking Jost's diversities
using Diversity.Jost

@test jostD == qD
@test jostβ == jostbeta
@test_approx_eq jostbeta(communities, 1) B̄(communities, 1)
@test_approx_eq jostbeta(allthesame, qs) ones(qs)

# Checking Hill numbers
using Diversity.Hill

@test hillnumber == qD
