module ValidateDistances
using Base.Test

using Diversity
using Distances
using Diversity.ShortNames

@testset "RenyiDivergence" begin
    @testset "Random RenyiDivergence $i" for i in 1:20
        types = rand(1:(i*10))
        sc = rand(1:(i*10))
        pops = rand(types, sc)
        # Make sure not to remove all of the non-zeros from any column
        for i in 1:sc
            pops[pops[:, i] .< median(pops[:, i])/2, i] = 0.0
        end
        pops /= sum(pops)
        meta = Metacommunity(pops)
        metapop = Diversity.vectorise(getmetaordinariness!(meta))
        weight_row = getweight(meta)
        weights = Diversity.vectorise(weight_row)
        pops_norm = pops ./ weight_row
        nb = β̄(meta)
        rb = ρ̄(meta)
        n = β(meta)
        r = ρ(meta)
        for q in [rand(7)*10..., 0, 1, Inf]
            @test exp.(colwise(RenyiDivergence(q),
                               pops_norm, metapop)) ≈ subdiv(nb, q)[:diversity]
            @test exp.(-colwise(RenyiDivergence(q),
                                pops_norm, metapop)) ≈ subdiv(rb, q)[:diversity]
            @test exp.(colwise(RenyiDivergence(q),
                               pops_norm, metapop)) ≈ subdiv(n, q)[:diversity] ./ weights
            @test exp.(-colwise(RenyiDivergence(q),
                                pops_norm, metapop)) ≈ subdiv(r, q)[:diversity] .* weights
        end
    end
end

end
