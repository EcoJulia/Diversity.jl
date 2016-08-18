using Diversity
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

include("EffectiveNumbers.jl")
include("Supercommunity.jl")
include("DiversityMeasure.jl")
include("Hill.jl")
include("Jost.jl")
include("Ecology.jl")

# Need to check the diversity() function
#@test diversity(Diversity.ρ̄, allthesame, qs, Z, false, false, true) ≈ colweights
#@test diversity((Diversity.ρ̄, :sub),
#                          communities, qs, Z) ≈ subcommunityrhobar(communities, qs, Z)
#@test diversity((Diversity.ρ̄, :super), communities, qs, Z) ≈ supercommunityRbar(communities, qs, Z)
#@test diversity((Diversity.ᾱ, Set([:super, :sub])), communities, qs, Z)[1] ≈ supercommunityAbar(communities, qs, Z)
#@test diversity((Diversity.α, Set([:super, :sub, :weights])), allthesame, qs, Z)[2] ≈ subcommunityalpha(allthesame, qs, Z)
#ed, cd, w = diversity((Diversity.γ, Set([:super, :sub, :weights])), communities, qs, Z)
#@test diversity((Diversity.γ̄, :super), communities, qs, Z) ≈ ed 
#@test diversity((Diversity.γ̄, Set([:sub, :weights])), communities, qs, Z)[1] ≈ cd
#@test diversity((Diversity.α, Set([:super, :sub, :weights])), allthesame, qs, Z)[3] ≈ colweights

#using Diversity.contributions
# Now some even communities, should see that raw and normalised
# diversities are the same
#smoothed = communities ./ mapslices(sum, communities, 1);
#smoothed /= numcommunities;
# Just for completeness, check one for q=-Inf - we currently have no use for this, but it is coded.
#@test contributions(Diversity.α, smoothed, [-Inf, 0, 1, 2, 3, 4, 5, Inf], true) ≈ contributions(Diversity.ᾱ, smoothed, [-Inf, 0, 1, 2, 3, 4, 5, Inf], true)
#@test contributions(Diversity.ρ, smoothed, qs, true) ≈ contributions(Diversity.ρ̄, smoothed, qs, true)
#@test contributions(Diversity.γ, smoothed, qs, true) ≈ contributions(Diversity.γ̄, smoothed, qs, true)
#@test contributions(Diversity.α, smoothed, qs, false) ≈ contributions(Diversity.ᾱ, smoothed, qs, false)
#@test contributions(Diversity.ρ, smoothed, qs, false) ≈ contributions(Diversity.ρ̄, smoothed, qs, false)
#@test contributions(Diversity.γ, smoothed, qs, false) ≈ contributions(Diversity.γ̄, smoothed, qs, false)

#@test contributions(Diversity.α, smoothed, qs, true) ≈ contributions(Diversity.α, smoothed, qs, false) * numcommunities
#@test contributions(Diversity.ρ, smoothed, qs, true) ≈ contributions(Diversity.ρ, smoothed, qs, false) * numcommunities
#@test contributions(Diversity.γ, smoothed, qs, true) ≈ contributions(Diversity.γ, smoothed, qs, false) * numcommunities
#@test contributions(Diversity.ᾱ, smoothed, qs, true) ≈ contributions(Diversity.ᾱ, smoothed, qs, false) * numcommunities
#@test contributions(Diversity.ρ̄, smoothed, qs, true) ≈ contributions(Diversity.ρ̄, smoothed, qs, false) * numcommunities
#@test contributions(Diversity.γ̄, smoothed, qs, true) ≈ contributions(Diversity.γ̄, smoothed, qs, false) * numcommunities
