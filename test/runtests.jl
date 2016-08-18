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
include("Supercommunity.jl")
include("DiversityMeasure.jl")
include("Hill.jl")
include("Jost.jl")
include("Ecology.jl")

# Sub-community alpha diversities
communities = rand(numspecies, numcommunities);
communities /= sum(communities);
@test_approx_eq subdiv(NormalisedAlpha(Ecosystem(communities)), 0) numspecies * ones((1, size(communities, 2)))
@test_approx_eq subdiv(NormalisedAlpha(Ecosystem(communities)),
                       [0])[1] numspecies * ones((1, size(communities, 2)))
sna = subdiv(NormalisedAlpha(Ecosystem(communities, Z1)),
             [0, 1, 2, Inf])
@test length(sna) == 4
for i in eachindex(sna)
    @test_approx_eq sna[i] ones(size(communities, 2))
end

@test_approx_eq subdiv(RawAlpha(Ecosystem(communities)), 0) numspecies * mapslices(v -> 1. / sum(v), communities, 1)

even = ones((numspecies, numcommunities)) / (numspecies * numcommunities);
qs = [0, 1, 2, 3, 4, 5, 6, Inf];
@test_approx_eq superdiv(NormalisedAlpha(Ecosystem(even)), qs) numspecies * ones((1, length(qs)))
@test_approx_eq superdiv(RawAlpha(Ecosystem(even)), qs) numspecies * numcommunities * ones((1, length(qs)))

probs = reshape(mapslices(sum, communities, 2), (size(communities, 1)));
@test_approx_eq superdiv(Gamma(Ecosystem(communities)), qs) qD(probs, qs)
@test_approx_eq superdiv(Gamma(Ecosystem(communities, Z1)), qs) qDZ(probs, qs, Z1)

Z = rand(numspecies, numspecies);
@test_approx_eq superdiv(Gamma(Ecosystem(communities, Z)), qs) qDZ(probs, qs, Z)

colweights = rand(numcommunities);
colweights /= sum(colweights);
allthesame = probs * colweights';
@test_approx_eq superdiv(RawBeta(Ecosystem(allthesame, Z)), qs) 1 ./ qD(colweights, 2 - qs)
@test_approx_eq superdiv(NormalisedBeta(Ecosystem(allthesame, Z)), qs) ones((1, length(qs)))
@test_approx_eq superdiv(NormalisedRho(Ecosystem(allthesame, Z)), qs) ones((1, length(qs)))
@test_approx_eq superdiv(RawRho(Ecosystem(allthesame, Z)), qs) qD(colweights, qs)

communitylist = rand(1:numcommunities, numspecies)
distinct = zeros(Float64, (numspecies, numcommunities))
for i in 1:numspecies
    distinct[i, communitylist[i]] = weights[i]
end

@test_approx_eq superdiv(RawRho(Ecosystem(distinct)), qs) ones((1, length(qs)))
subnr = subdiv(NormalisedRho(Ecosystem(distinct)), qs)
for i in eachindex(qs)
    @test_approx_eq subnr[i] sum(distinct, 1)
end
@test_approx_eq superdiv(NormalisedBeta(Ecosystem(distinct)), qs) qD(reshape(sum(distinct, 1), numcommunities), qs)
@test_approx_eq superdiv(RawBeta(Ecosystem(distinct)), qs) ones((1, length(qs)))

# Need to check the diversity() function
#@test_approx_eq diversity(Diversity.ρ̄, allthesame, qs, Z, false, false, true) colweights
#@test_approx_eq diversity((Diversity.ρ̄, :sub),
#                          communities, qs, Z) subcommunityrhobar(communities, qs, Z)
#@test_approx_eq diversity((Diversity.ρ̄, :super), communities, qs, Z) supercommunityRbar(communities, qs, Z)
#@test_approx_eq diversity((Diversity.ᾱ, Set([:super, :sub])), communities, qs, Z)[1] supercommunityAbar(communities, qs, Z)
#@test_approx_eq diversity((Diversity.α, Set([:super, :sub, :weights])), allthesame, qs, Z)[2] subcommunityalpha(allthesame, qs, Z)
#ed, cd, w = diversity((Diversity.γ, Set([:super, :sub, :weights])), communities, qs, Z)
#@test_approx_eq diversity((Diversity.γ̄, :super), communities, qs, Z) ed 
#@test_approx_eq diversity((Diversity.γ̄, Set([:sub, :weights])), communities, qs, Z)[1] cd
#@test_approx_eq diversity((Diversity.α, Set([:super, :sub, :weights])), allthesame, qs, Z)[3] colweights

#using Diversity.contributions
# Now some even communities, should see that raw and normalised
# diversities are the same
smoothed = communities ./ mapslices(sum, communities, 1);
smoothed /= numcommunities;
# Just for completeness, check one for q=-Inf - we currently have no use for this, but it is coded.
#@test_approx_eq contributions(Diversity.α, smoothed, [-Inf, 0, 1, 2, 3, 4, 5, Inf], true) contributions(Diversity.ᾱ, smoothed, [-Inf, 0, 1, 2, 3, 4, 5, Inf], true)
#@test_approx_eq contributions(Diversity.ρ, smoothed, qs, true) contributions(Diversity.ρ̄, smoothed, qs, true)
#@test_approx_eq contributions(Diversity.γ, smoothed, qs, true) contributions(Diversity.γ̄, smoothed, qs, true)
#@test_approx_eq contributions(Diversity.α, smoothed, qs, false) contributions(Diversity.ᾱ, smoothed, qs, false)
#@test_approx_eq contributions(Diversity.ρ, smoothed, qs, false) contributions(Diversity.ρ̄, smoothed, qs, false)
#@test_approx_eq contributions(Diversity.γ, smoothed, qs, false) contributions(Diversity.γ̄, smoothed, qs, false)

#@test_approx_eq contributions(Diversity.α, smoothed, qs, true) contributions(Diversity.α, smoothed, qs, false) * numcommunities
#@test_approx_eq contributions(Diversity.ρ, smoothed, qs, true) contributions(Diversity.ρ, smoothed, qs, false) * numcommunities
#@test_approx_eq contributions(Diversity.γ, smoothed, qs, true) contributions(Diversity.γ, smoothed, qs, false) * numcommunities
#@test_approx_eq contributions(Diversity.ᾱ, smoothed, qs, true) contributions(Diversity.ᾱ, smoothed, qs, false) * numcommunities
#@test_approx_eq contributions(Diversity.ρ̄, smoothed, qs, true) contributions(Diversity.ρ̄, smoothed, qs, false) * numcommunities
#@test_approx_eq contributions(Diversity.γ̄, smoothed, qs, true) contributions(Diversity.γ̄, smoothed, qs, false) * numcommunities
