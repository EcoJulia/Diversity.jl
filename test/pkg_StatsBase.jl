module ValidateEntropies
using Compat.Test
using Compat
using Compat.Statistics

using Diversity
using Diversity.ShortNames
using StatsBase

@testset "renyientropy" begin
    @testset "Random renyientropy $i" for i in 1:20
        types = rand(1:(i*10))
        pop = rand(types)
        pop[pop .< Compat.Statistics.median(pop)/2] .= 0.0
        pop ./= sum(pop)
        meta = Metacommunity(pop)
        g = Γ(meta)
        na = ᾱ(meta)
        ra = α(meta)
        qs = [rand(7)*10..., 0, 1, Inf]
        entropies = map(q -> renyientropy(pop, q), qs)
        @test log.(subdiv(na, qs)[:diversity]) ≈ entropies
        @test log.(subdiv(ra, qs)[:diversity]) ≈ entropies
        @test log.(subdiv(g, qs)[:diversity]) ≈ entropies
        @test log.(qD(pop, qs)) ≈ entropies
        @test log.(qDZ(pop, qs, UniqueTypes(types))) ≈ entropies
    end
end

end
