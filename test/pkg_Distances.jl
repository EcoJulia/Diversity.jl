module ValidateDistances
using Compat.Test
using Compat
using Compat.Statistics

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
            vals = pops[:, j] .< Compat.Statistics.median(pops[:, j])/2
            for k in 1:types
                if vals[k]
                    pops[k, j] = 0.0
                end
            end
        end
        pops /= sum(pops)
        meta = Metacommunity(pops)
        metapop = getmetaordinariness!(meta)
        weights = getweight(meta)
        pops_norm = pops ./ weights'
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
