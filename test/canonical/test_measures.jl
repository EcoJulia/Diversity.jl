# SPDX-License-Identifier: BSD-2-Clause
#
# Canonical results for the **measures themselves**: a fixed abundance matrix and a fixed similarity
# matrix in, every measure at every scale out. Nothing is read and nothing is random, so every number
# here is a pure function of the arithmetic in `src/` - which makes this the file that isolates a
# change in a *measure* from a change in how a *type* builds its similarity matrix (`test_types.jl`).

module CanonicalMeasures

using Test
using Diversity
using Diversity.ShortNames

include("canonical.jl")
using .Canonical

# A deliberately lopsided metacommunity: unequal subcommunity weights, and a zero in each of two
# subcommunities, so the measures are not accidentally symmetric in a way that would hide a change.
const POP = [2.0 1.0 0.0
             1.0 3.0 1.0
             0.0 1.0 4.0]

# **Asymmetric on purpose.** Similarity need not be symmetric - `run_rcall.jl` cross-validates
# both cases against rdiversity - so blessing a symmetric matrix here would leave the asymmetric path
# unpinned, which is the one where a transposed index goes unnoticed.
const ZASYM = [1.0 0.5 0.0
               0.2 1.0 0.4
               0.0 0.3 1.0]

const QS = [0, 1, 2, Inf]
qlabel(q) = isinf(q) ? "qInf" : "q$(Int(q))"

const MEASURES = ["RawAlpha" => RawAlpha, "NormalisedAlpha" => NormalisedAlpha,
    "RawBeta" => RawBeta, "NormalisedBeta" => NormalisedBeta,
    "RawRho" => RawRho, "NormalisedRho" => NormalisedRho,
    "Gamma" => Gamma]

@testset "canonical: measures" begin
    @testset "$tname types" for (tname, types) in ["unique" => UniqueTypes(3),
                                    "asymmetricz" => GeneralTypes(ZASYM)]
        meta = Metacommunity(POP, types)

        # Both scales are blessed, not just the metacommunity. They aggregate with *different*
        # power-mean orders - opposite ones for the relative-entropy measures - so a subcommunity
        # vector can move while its metacommunity summary does not.
        for (mname, measure) in MEASURES
            dm = measure(meta)
            for q in QS
                blessed("measures/$tname/$mname/sub_$(qlabel(q))",
                        subdiv(dm, q)[!, :diversity])
                blessed("measures/$tname/$mname/meta_$(qlabel(q))",
                        metadiv(dm, q)[1, :diversity])
            end
        end

        # The properties that must hold whatever the blessed numbers are. Re-blessing silences the
        # values above; it must never be able to silence these.
        w = getweight(meta)
        for q in QS
            # β̄ and ρ̄ are exact reciprocals individually, and the opposite aggregation orders carry
            # that through to the subcommunity - but *not* to the metacommunity, where both use
            # order 1 - q. Asserting the second would be wrong; asserting the first pins the design.
            @test subdiv(β̄(meta), q)[!, :diversity] .*
                  subdiv(ρ̄(meta), q)[!, :diversity] ≈ ones(3)
            # raw vs normalised differ by exactly the subcommunity weight
            @test subdiv(α(meta), q)[!, :diversity] ≈
                  subdiv(ᾱ(meta), q)[!, :diversity] ./ w
            @test all(isfinite, metadiv(Γ(meta), q)[!, :diversity])
            @test all(>(0), subdiv(Γ(meta), q)[!, :diversity])
        end
    end

    # Note: Numbers equivalence: n equally-abundant, wholly distinct types must give exactly n at every
    # q. It is the defining property of the whole framework and the cheapest thing in this file.
    @testset "numbers equivalence" begin
        even = Metacommunity(fill(1 / 5, 5))
        for q in QS
            @test metadiv(Γ(even), q)[1, :diversity] ≈ 5
            @test metadiv(ᾱ(even), q)[1, :diversity] ≈ 5
        end
    end
end

end
