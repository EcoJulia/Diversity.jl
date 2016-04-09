using Diversity
using Base.Test

numbers = [1., 2, 4, 8, 16];
numspecies = 100;
fragments = rand(numspecies);
weights = rand(numspecies);
weights /= sum(weights);
Z1 = ones(typeof(weights[1]), (length(weights), length(weights)));
numcommunities = 8;
manyweights = rand(numspecies, numcommunities);
manyweights *= diagm(reshape(mapslices(v -> 1. / sum(v), manyweights, 1),
                             (numcommunities)));

include("EffectiveNumbers.jl")

# Sub-community alpha diversities
communities = rand(numspecies, numcommunities);
communities /= sum(communities);
@test_approx_eq Dᾱ(communities, 0) numspecies * ones((1, size(communities, 2)))
@test_approx_eq subcommunityalphabar(communities,
                                     [0]) numspecies * ones((1, size(communities, 2)))
@test_approx_eq Dᾱ(communities,
                  [0, 1, 2, Inf], Z1) ones((4, size(communities, 2)))

@test_approx_eq Dα(communities, 0) numspecies * mapslices(v -> 1. / sum(v),
                                                  communities, 1)

even = ones((numspecies, numcommunities)) / (numspecies * numcommunities);
qs = [0, 1, 2, 3, 4, 5, 6, Inf];
@test_approx_eq DĀ(even, qs) numspecies * ones((1, length(qs)))
@test_approx_eq DA(even, qs) numspecies * numcommunities * ones((1, length(qs)))

probs = reshape(mapslices(sum, communities, 2), (size(communities, 1)));
@test_approx_eq DG(communities, qs) DḠ(communities, qs)
@test_approx_eq DG(communities, qs) qD(probs, qs)
@test_approx_eq DG(communities, qs, Z1) qDZ(probs, qs, Z1)

Z = rand(numspecies, numspecies);
@test_approx_eq DG(communities, qs, Z) qDZ(probs, qs, Z)

colweights = rand(numcommunities);
colweights /= sum(colweights);
allthesame = probs * colweights';
@test_approx_eq DB(allthesame, qs, Z) 1 ./ qD(colweights, 2 - qs)
@test_approx_eq DB̄(allthesame, qs, Z) ones((1, length(qs)))
@test_approx_eq DR̄(allthesame, qs, Z) ones((1, length(qs)))
@test_approx_eq DR(allthesame, qs) qD(colweights, qs)

communitylist = rand(1:numcommunities, numspecies)
distinct = zeros(Float64, (numspecies, numcommunities))
for i in 1:numspecies
    distinct[i, communitylist[i]] = weights[i]
end

@test_approx_eq DR(distinct, qs) ones((1, length(qs)))
@test_approx_eq Dρ̄(distinct, qs) repeat(sum(distinct, 1), inner = [length(qs), 1])
@test_approx_eq DB̄(distinct, qs) qD(reshape(sum(distinct, 1), numcommunities), qs)
@test_approx_eq DB(distinct, qs) ones((1, length(qs)))

# Need to check the diversity() function
@test_approx_eq diversity(DR̄, allthesame, qs, Z, false, false, true) colweights
@test_approx_eq diversity(Dρ̄, communities, qs, Z,
                          false, true, false) Dρ̄(communities, qs, Z)
@test_approx_eq diversity(Dρ̄, communities, qs, Z,
                          true, false, false) DR̄(communities, qs, Z)
@test_approx_eq diversity(Dᾱ, communities, qs, Z,
                          true, true, false)[1] DĀ(communities, qs, Z)
@test_approx_eq diversity(Dα, allthesame, qs, Z,
                          true, true, true)[2] Dα(allthesame, qs, Z)
ed, cd, w = diversity(Dγ, communities, qs, Z, true, true, true)
@test_approx_eq diversity(Dγ̄, communities, qs, Z, true, false, true)[1] ed 
@test_approx_eq diversity(Dγ̄, communities, qs, Z, false, true, true)[1] cd
@test_approx_eq diversity(Dα, allthesame, qs, Z, true, true, true)[3] colweights

using Diversity.contributions
# Now some even communities, should see that raw and normalised
# diversities are the same
smoothed = communities ./ mapslices(sum, communities, 1);
smoothed /= numcommunities;
# Just for completeness, check one for q=-Inf - we currently have no use for this, but it is coded.
@test_approx_eq contributions(Dα, smoothed, [-Inf, 0, 1, 2, 3, 4, 5, Inf], true) contributions(Dᾱ, smoothed, [-Inf, 0, 1, 2, 3, 4, 5, Inf], true)
@test_approx_eq contributions(Dρ, smoothed, qs, true) contributions(Dρ̄, smoothed, qs, true)
@test_approx_eq contributions(Dγ, smoothed, qs, true) contributions(Dγ̄, smoothed, qs, true)
@test_approx_eq contributions(Dα, smoothed, qs, false) contributions(Dᾱ, smoothed, qs, false)
@test_approx_eq contributions(Dρ, smoothed, qs, false) contributions(Dρ̄, smoothed, qs, false)
@test_approx_eq contributions(Dγ, smoothed, qs, false) contributions(Dγ̄, smoothed, qs, false)

@test_approx_eq contributions(Dα, smoothed, qs, true) contributions(Dα, smoothed, qs, false) * numcommunities
@test_approx_eq contributions(Dρ, smoothed, qs, true) contributions(Dρ, smoothed, qs, false) * numcommunities
@test_approx_eq contributions(Dγ, smoothed, qs, true) contributions(Dγ, smoothed, qs, false) * numcommunities
@test_approx_eq contributions(Dᾱ, smoothed, qs, true) contributions(Dᾱ, smoothed, qs, false) * numcommunities
@test_approx_eq contributions(Dρ̄, smoothed, qs, true) contributions(Dρ̄, smoothed, qs, false) * numcommunities
@test_approx_eq contributions(Dγ̄, smoothed, qs, true) contributions(Dγ̄, smoothed, qs, false) * numcommunities

# Looking at relations to historical measures, updated with similarity
# and partitioning
using Diversity.Ecology

@test_approx_eq richness(communities) Dᾱ(communities, 0)
@test_approx_eq generalisedrichness(DḠ, communities, Z1) 1

@test_approx_eq shannon(communities) log(Dᾱ(communities, 1))
@test_approx_eq generalisedshannon(DḠ, communities, Z1) 0

@test_approx_eq simpson(communities) Dᾱ(communities, 2) .^ -1
@test_approx_eq generalisedsimpson(DḠ, communities, Z1) 1

@test_approx_eq jaccard([1 0 0; 0 1 1.]') 0
@test_approx_eq jaccard([1 0 1; 0 1 1.]') 1 / 3
@test_throws ErrorException jaccard([1 1 0; 0 1 1; 1 1 1.])

@test_approx_eq generalisedjaccard([1 0 1; 0 1 1.]', [0, Inf]) [1/3, 1]
@test_approx_eq generalisedjaccard([1 1 1; 1 1 1.]', [0, 1]) [1, 1]

# Checking Jost's diversities
using Diversity.Jost

@test jostβ == jostbeta
@test_approx_eq jostbeta(communities, 1) 1 ./ DR̄(communities, 1)
@test_approx_eq jostbeta(allthesame, qs) ones(qs)

## Check Jost's alpha diversity works for all the same subcommunity
@test jostα == jostalpha
@test_approx_eq jostalpha(allthesame, qs) DĀ(allthesame, qs)
## And for all different subcommunities and any subcommunities with the same sizes
evendistinct = mapslices((x) -> x / (sum(x) * numcommunities), distinct, 1)
@test_approx_eq jostalpha(evendistinct, qs) DĀ(evendistinct, qs)
@test_approx_eq jostalpha(smoothed, qs) DĀ(smoothed, qs)

# Checking Hill numbers
using Diversity.Hill

@test_approx_eq hillnumber(manyweights, 0) numspecies * ones((1, size(manyweights, 2)))
