# SPDX-License-Identifier: BSD-2-Clause

module TestDiversity
using Test

using Diversity
using Diversity.ShortNames

@testset "Package structure" begin
    # `path` locates files inside the package, defaulting to `test/`. It is what the extension tests
    # and the genetics documentation use to reach `test/data`.
    @test isfile(Diversity.path("runtests.jl"))
    @test isfile(Diversity.path("data", "biallelic.vcf"))
    @test isdir(Diversity.path("Diversity.jl"; dir = "src") |> dirname)

    # The submodules are part of the public structure, not implementation detail.
    for sub in (:API, :ShortNames, :Ecology, :Jost, :Hill)
        @test isdefined(Diversity, sub)
    end
end

@testset "ShortNames" begin
    # The unicode names are aliases, not separate implementations — if these ever diverge, the two
    # spellings of the same measure would silently disagree.
    @test α ≡ RawAlpha
    @test ᾱ ≡ NormalisedAlpha
    @test β ≡ RawBeta
    @test β̄ ≡ NormalisedBeta
    @test ρ ≡ RawRho
    @test ρ̄ ≡ NormalisedRho
    @test Γ ≡ Gamma

    # ⚠️ `γ` cannot be exported (Julia always resolves it as `ShortNames.γ`), so `Γ` is exported in
    # its place. That is deliberate, and this is the assertion that says so.
    @test Diversity.ShortNames.γ ≡ Gamma
    @test :Γ ∈ names(Diversity.ShortNames)
    @test :γ ∉ names(Diversity.ShortNames)
end

@testset "_dist2sim" begin
    dist = [0.0 4.0 3.0
            4.0 0.0 3.0
            3.0 3.0 0.0]

    # Linear with normalisation is the rdiversity default: Z = 1 - d / max(d).
    linear = Diversity._dist2sim(dist; transform = :linear, k = 1,
                                 normalise = true, max_d = maximum(dist))
    @test linear ≈ 1.0 .- dist ./ 4.0
    @test linear isa Matrix{Float64}

    # Exponential is the other transform rdiversity's dist2sim offers.
    expo = Diversity._dist2sim(dist; transform = :exponential, k = 1,
                               normalise = true, max_d = maximum(dist))
    @test expo ≈ exp.(-dist ./ 4.0)

    # Without normalisation the raw distances are used directly.
    raw = Diversity._dist2sim(dist; transform = :exponential, k = 1,
                              normalise = false, max_d = maximum(dist))
    @test raw ≈ exp.(-dist)

    # k scales the distance before the transform.
    @test Diversity._dist2sim(dist; transform = :exponential, k = 2,
                              normalise = false, max_d = maximum(dist)) ≈
          exp.(-2 .* dist)

    # ⚠️ The linear transform clamps at zero rather than going negative — a similarity below 0 would
    # be rejected by `GeneralTypes` and is meaningless anyway.
    clamped = Diversity._dist2sim(dist; transform = :linear, k = 10,
                                  normalise = true, max_d = maximum(dist))
    @test all(≥(0), clamped)
    @test clamped[1, 2] == 0.0
    @test all(≈(1.0), clamped[i, i] for i in 1:3)

    # An all-zero distance matrix has no maximum to normalise by; dividing anyway would give NaN.
    zeroes = zeros(2, 2)
    @test Diversity._dist2sim(zeroes; transform = :linear, k = 1,
                              normalise = true, max_d = 0.0) == ones(2, 2)

    # 🔴 The error branch. Reachable in normal use only by passing a bad `transform` to
    # `GeneticType`, so nothing else in the suite covers it.
    @test_throws ArgumentError Diversity._dist2sim(dist; transform = :quadratic,
                                                   k = 1, normalise = true,
                                                   max_d = maximum(dist))
end

@testset "Genetic stubs" begin
    # Both are declared method-less in the parent so the extensions can add the sole method — the
    # "extensions add, never overwrite" rule. Without an extension loaded they exist but do nothing.
    @test Diversity.GeneticType isa Function
    @test Diversity.vcf_dataframe isa Function
    @test Diversity.AbstractGenetic <: Diversity.API.AbstractTypes
end

end
