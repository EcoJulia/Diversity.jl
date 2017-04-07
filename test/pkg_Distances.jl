module ValidateDistances
using Base.Test

using Diversity
using Distances
using Diversity.ShortNames

@testset "RenyiDivergence" for i in 1:10
    types = rand(1:100)
    sc = rand(1:100)
    pops = rand(types, sc)
    # Make sure not to remove all of the non-zeros from any column
    for i in 1:sc
        pops[pops[:, i] .< median(pops[:, i])/2, i] = 0.0
    end
    pops /= sum(pops)
    meta = Metacommunity(pops)
    metapop = Diversity.vectorise(getmetaordinariness!(meta))
    weights = Diversity.vectorise(getweight(meta))
    nb = β̄(meta)
    rb = ρ̄(meta)
    n = β(meta)
    r = ρ(meta)
    for q in [rand(10)*10..., 0, 1, Inf]
        @test exp.(colwise(RenyiDivergence(q), pops, metapop)) ≈ subdiv(nb, q)[:diversity]
        @test exp.(-colwise(RenyiDivergence(q), pops, metapop)) ≈ subdiv(rb, q)[:diversity]
        @test exp.(colwise(RenyiDivergence(q), pops, metapop)) ≈ subdiv(n, q)[:diversity] ./ weights
        @test exp.(-colwise(RenyiDivergence(q), pops, metapop)) ≈ subdiv(r, q)[:diversity] .* weights
    end
end

end
