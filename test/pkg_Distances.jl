module ValidateDistances
using Base.Test

using Diversity
using Diversity.ShortNames
using Distances

@testset "RenyiDivergence" begin
    @testset "Random RenyiDivergence $i" for i in 1:20
        types = rand(1:(i*10))
        sc = rand(1:(i*10))
        pops = rand(types, sc)
        # Make sure not to remove all of the non-zeros from any column
        for j in 1:sc
            pops[pops[:, j] .< median(pops[:, j])/2, j] = 0.0
        end
        pops /= sum(pops)
        meta = Metacommunity(pops)
        metapop = Diversity.vectorise(getmetaordinariness!(meta))
        weight_row = getweight(meta)
        weights = Diversity.vectorise(weight_row)
        pops_norm = pops ./ weight_row
        nb = β̄(meta)
        nr = ρ̄(meta)
        b = β(meta)
        r = ρ(meta)
        for q in [rand(7)*10..., 0, 1, Inf]
            @test exp.(colwise(RenyiDivergence(q),
                               pops_norm, metapop)) ≈ subdiv(nb, q)[:diversity]
            @test exp.(-colwise(RenyiDivergence(q),
                                pops_norm, metapop)) ≈ subdiv(nr, q)[:diversity]
            @test exp.(colwise(RenyiDivergence(q),
                               pops_norm, metapop)) ≈ subdiv(b, q)[:diversity] ./ weights
            @test exp.(-colwise(RenyiDivergence(q),
                                pops_norm, metapop)) ≈ subdiv(r, q)[:diversity] .* weights
        end
    end
end

end
