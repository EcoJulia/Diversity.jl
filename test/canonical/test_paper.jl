# SPDX-License-Identifier: BSD-2-Clause
#
# The package checked against the **published** values in Reeve et al, "How to partition diversity"
# (arXiv:1404.6520), whose appendices work several communities right through and state every measure.
#
# 🔴 **Nothing here is blessed, and nothing here may become blessed.** Every other file in this
# directory records what the package produced, and re-blessing rewrites it. These are plain `@test`s
# against numbers the paper states, so no re-blessing can touch them: if one fails, either the
# package has stopped implementing the specification or the specification has changed, and both are
# findings rather than something to record and move past.
#
# ⚠️ Examples are named rather than numbered — supplementary section numbers move between versions of
# a paper, and a pointer a reader cannot follow is worse than none.

module CanonicalPaper

using Test
using Diversity
using Diversity.ShortNames

const QS = [0, 1, 2, Inf]

@testset "paper: the sexual contact network" begin
    # Six individuals as types, similarity being the chance of transmission, divided into three
    # males and three females. Transmission is likelier between subcommunities than within them, so
    # every type is *less* similar to its own subcommunity than to the other.
    Z = [1.0 0.5 0.5 0.7 0.7 0.7
         0.5 1.0 0.5 0.7 0.7 0.7
         0.5 0.5 1.0 0.7 0.7 0.7
         0.7 0.7 0.7 1.0 0.5 0.5
         0.7 0.7 0.7 0.5 1.0 0.5
         0.7 0.7 0.7 0.5 0.5 1.0]
    P = [1 0; 1 0; 1 0; 0 1; 0 1; 0 1] ./ 6
    meta = Metacommunity(P, Z)

    # ⭐ Every ordinariness here is 4.1/6, which is what makes every measure independent of q and the
    # published values exact rather than rounded: 2/4.1, 4/4.1 and 6/4.1 in place of the 0.488, 0.976
    # and 1.463 the paper prints to three decimal places.
    @testset "q = $q" for q in [0.5, 1, 2, 5]
        @test raw_sub_alpha(meta, q)[!, :diversity] ≈ [3.0, 3.0]
        @test norm_sub_alpha(meta, q)[!, :diversity] ≈ [1.5, 1.5]
        @test raw_sub_rho(meta, q)[!, :diversity] ≈ [2.05, 2.05]
        @test raw_sub_beta(meta, q)[!, :diversity] ≈ [2 / 4.1, 2 / 4.1]
        @test norm_sub_rho(meta, q)[!, :diversity] ≈ [1.025, 1.025]
        @test norm_sub_beta(meta, q)[!, :diversity] ≈ [4 / 4.1, 4 / 4.1]
        @test sub_gamma(meta, q)[!, :diversity] ≈ [6 / 4.1, 6 / 4.1]

        # The metacommunity measures take the same values as the subcommunity ones here.
        @test raw_meta_alpha(meta, q)[1, :diversity] ≈ 3.0
        @test norm_meta_alpha(meta, q)[1, :diversity] ≈ 1.5
        @test raw_meta_rho(meta, q)[1, :diversity] ≈ 2.05
        @test raw_meta_beta(meta, q)[1, :diversity] ≈ 2 / 4.1
        @test norm_meta_rho(meta, q)[1, :diversity] ≈ 1.025
        @test norm_meta_beta(meta, q)[1, :diversity] ≈ 4 / 4.1
        @test meta_gamma(meta, q)[1, :diversity] ≈ 6 / 4.1
    end

    # ⭐ This is why the example is in the paper at all: with a general similarity matrix the bounds
    # that hold when types are wholly distinct all fail. Representativeness exceeds 1 and
    # distinctiveness falls below it, neither of which can happen when Z = I.
    @test all(>(1), norm_sub_rho(meta, 1)[!, :diversity])
    @test all(<(1), norm_sub_beta(meta, 1)[!, :diversity])

    # ⚠️ And it is not a pathological matrix — it satisfies the triangle inequality — so no weaker
    # assumption than Z = I rescues those bounds.
    @test all(Z[i, j] * Z[j, k] ≤ Z[i, k] + eps()
              for i in 1:6, j in 1:6, k in 1:6)
end

@testset "paper: representativeness is a proportion" begin
    # All types equally abundant, each subcommunity holding an equal share of them: a subcommunity
    # holding a fraction r of the types has representativeness exactly r.
    P = [1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1] ./ 6
    meta = Metacommunity(P)
    for q in QS
        @test norm_sub_rho(meta, q)[!, :diversity] ≈ fill(1 / 3, 3)
        @test norm_meta_rho(meta, q)[1, :diversity] ≈ 1 / 3
    end
end

@testset "paper: subcommunity size does not change gamma" begin
    # Two subcommunities of types that are equally abundant in the metacommunity, the first holding
    # twice as many of them: the first is twice as diverse in isolation, but each individual
    # contributes the same to the metacommunity.
    P = [1 0; 1 0; 1 0; 1 0; 0 1; 0 1] ./ 6
    meta = Metacommunity(P)
    for q in QS
        alphas = norm_sub_alpha(meta, q)[!, :diversity]
        gammas = sub_gamma(meta, q)[!, :diversity]
        @test alphas ≈ [4.0, 2.0]
        @test alphas[1] ≈ 2 * alphas[2]
        @test gammas ≈ [6.0, 6.0]
        @test gammas[1] ≈ gammas[2]
    end
end

@testset "paper: rarity does change gamma" begin
    # Two equally sized subcommunities, equally diverse in isolation, but the first one's types are
    # three times rarer in the metacommunity — so it contributes three times as much per individual.
    P = [1 0 0; 1 0 0; 0 1 2; 0 1 2] ./ 8
    meta = Metacommunity(P)
    @test getweight(meta)[1] ≈ getweight(meta)[2]
    for q in QS
        alphas = norm_sub_alpha(meta, q)[!, :diversity]
        gammas = sub_gamma(meta, q)[!, :diversity]
        @test alphas[1] ≈ alphas[2]
        @test gammas[1] ≈ 3 * gammas[2]
    end
end

@testset "paper: a metacommunity that is a single subcommunity" begin
    # With one subcommunity every beta measure is 1 and the alphas and gamma all collapse onto the
    # diversity of the community itself.
    p = [0.5, 0.3, 0.2]
    meta = Metacommunity(p, UniqueTypes(3), Onecommunity())
    for q in QS
        div = meta_gamma(meta, q)[1, :diversity]
        @test raw_sub_alpha(meta, q)[1, :diversity] ≈ div
        @test norm_sub_alpha(meta, q)[1, :diversity] ≈ div
        @test sub_gamma(meta, q)[1, :diversity] ≈ div
        for measure in (raw_sub_rho, raw_sub_beta, norm_sub_rho, norm_sub_beta)
            @test measure(meta, q)[1, :diversity] ≈ 1.0
        end
    end
end

@testset "paper: a well-mixed metacommunity" begin
    # Every subcommunity with the metacommunity's own composition. Then the normalised measures do
    # not see the division at all: representativeness and the effective number of subcommunities are
    # both 1, and average subcommunity diversity equals metacommunity diversity.
    p = [0.5, 0.3, 0.2]
    w = [0.6, 0.4]
    meta = Metacommunity(p * w')
    for q in QS
        @test norm_sub_rho(meta, q)[!, :diversity] ≈ ones(2)
        @test norm_sub_beta(meta, q)[!, :diversity] ≈ ones(2)
        @test norm_meta_rho(meta, q)[1, :diversity] ≈ 1.0
        @test norm_meta_beta(meta, q)[1, :diversity] ≈ 1.0
        @test norm_meta_alpha(meta, q)[1, :diversity] ≈
              meta_gamma(meta, q)[1, :diversity]
        # Raw redundancy is the Hill number of the subcommunity weights.
        @test raw_meta_rho(meta, q)[1, :diversity] ≈ qD(w, q)
    end
end

@testset "paper: a naive-community metacommunity" begin
    # Subcommunities with no shared types and no similarity between them. Distinctiveness is then at
    # its maximum of 1, the effective number of subcommunities is the Hill number of their sizes,
    # and metacommunity gamma equals raw metacommunity alpha.
    P = [3 0; 1 0; 0 2; 0 2] ./ 8
    meta = Metacommunity(P)
    w = getweight(meta)
    for q in QS
        @test raw_sub_beta(meta, q)[!, :diversity] ≈ ones(2)
        @test raw_sub_rho(meta, q)[!, :diversity] ≈ ones(2)
        @test raw_meta_beta(meta, q)[1, :diversity] ≈ 1.0
        @test raw_meta_rho(meta, q)[1, :diversity] ≈ 1.0
        @test norm_meta_beta(meta, q)[1, :diversity] ≈ qD(w, q)
        @test meta_gamma(meta, q)[1, :diversity] ≈
              raw_meta_alpha(meta, q)[1, :diversity]
        # Normalised subcommunity beta is the reciprocal of the subcommunity's size.
        @test norm_sub_beta(meta, q)[!, :diversity] ≈ 1 ./ w
    end
end

end
