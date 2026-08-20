# SPDX-License-Identifier: BSD-2-Clause

module TestDiversityMeasure
using Test
using LinearAlgebra

using Diversity
using Diversity.ShortNames
using DataFrames
using Tables
using EcoBase: getcoords
using Plots

pop = [3, 3, 4]
pop = pop / sum(pop)
oc = Onecommunity()
sim = [1.0 0 0; 1.0 1.0 0.0; 1.0 1.0 1.0]
ms = GeneralTypes(sim)
meta = Metacommunity(pop, ms, oc)
meta1 = Metacommunity(pop)

ab3 = [1.0 2; 3 0; 0 4]'
ab3 = ab3 / sum(ab3)
sc = Subcommunities(3)
sp = Species(2)
meta2 = Metacommunity(ab3, sp, sc)
nab = NormalisedAlpha(meta2)

@testset "Diversity measures" begin
    diversities = [RawAlpha, NormalisedAlpha, RawBeta, NormalisedBeta,
        RawRho, NormalisedRho, Gamma]
    shortds = [α, ᾱ, β, β̄, ρ, ρ̄, Γ]
    chars = ["α", "ᾱ", "β", "β̄", "ρ", "ρ̄", "γ"]
    asciis = ["RawAlpha", "NormalisedAlpha",
        "RawBeta", "NormalisedBeta",
        "RawRho", "NormalisedRho", "Gamma"]
    # These are the paper's own descriptions of the measures, at subcommunity level — which is
    # the level `getFullName`'s only consumer, the plot recipe, works at.
    fulls = ["estimate of naive-community metacommunity diversity",
        "diversity of subcommunity in isolation",
        "distinctiveness",
        "estimate of effective number of distinct subcommunities",
        "redundancy", "representativeness",
        "contribution per individual toward metacommunity diversity"]
    for i in axes(diversities, 1)
        @test diversities[i] == shortds[i]
        @test getName(diversities[i](meta)) == chars[i]
        @test getASCIIName(diversities[i](meta2)) == asciis[i]
        @test getFullName(diversities[i](meta1)) == fulls[i]
    end

    # The descriptive aliases are the names the papers use, and are the same types — so a measure
    # reached through one spelling must be identical to the same measure reached through another.
    @test Diversity.Distinctiveness ≡ RawBeta
    @test Diversity.Redundancy ≡ RawRho
    @test Diversity.Representativeness ≡ NormalisedRho
    @test getFullName(Diversity.Representativeness(meta)) ==
          "representativeness"

    # `getASCIIName` strips the module prefix and the type parameters, so it names the *measure*
    # rather than the concrete parameterisation — which is what the output DataFrame carries.
    @test !occursin("Diversity.", getASCIIName(Gamma(meta)))
    @test !occursin("{", getASCIIName(Gamma(meta)))
end

numbers = [1.0, 2, 4, 8, 16]
numspecies = 100
fragments = rand(numspecies)
weights = rand(numspecies)
weights /= sum(weights)
Z1 = ones(typeof(weights[1]), (length(weights), length(weights)))
numcommunities = 8
manyweights = rand(numspecies, numcommunities)
manyweights *= Diagonal(reshape(mapslices(v -> 1.0 / sum(v), manyweights;
                                          dims = 1), numcommunities))

@testset "inddiv / subdiv / metadiv" begin
    @test individualDiversity(nab, 0)[!, :diversity] ≈
          inddiv(nab, 0)[!, :diversity]
    @test individualDiversity(nab)(1)[!, :diversity] ≈
          inddiv(nab, 1)[!, :diversity]
    idnab = inddiv(nab, [2, 3])
    @test idnab[isapprox.(collect(idnab[!, :q]), 2), :diversity] ≈
          inddiv(nab, 2)[!, :diversity]
    allid = inddiv(meta2, Inf)
    @test allid[allid[!, :measure] .== "RawAlpha", :diversity] ≈
          inddiv(RawAlpha(meta2), Inf)[!, :diversity]

    @test subcommunityDiversity(nab, Inf)[!, :diversity] ≈
          subdiv(nab, Inf)[!, :diversity]
    sdnab = subdiv(nab, [4, 5])
    @test sdnab[isapprox.(collect(sdnab[!, :q]), 4), :diversity] ≈
          subdiv(nab, 4)[!, :diversity]
    allsd = subdiv(meta2, Inf)
    @test allsd[allsd[!, :measure] .== "NormalisedAlpha", :diversity] ≈
          subdiv(nab, Inf)[!, :diversity]
    allmd = subdiv(meta1, Inf)
    @test allmd[allmd[!, :measure] .== "NormalisedAlpha", :diversity] ≈
          metadiv(ᾱ(meta1), Inf)[!, :diversity]

    scg = subcommunityDiversity(Gamma(meta2))
    @test scg(1)[!, :diversity] ≈ scg(1.0)[!, :diversity]

    communities = rand(numspecies, numcommunities)
    communities /= sum(communities)
    @test subdiv(NormalisedAlpha(Metacommunity(communities)), 0)[!,
                                                                 :diversity] ≈
          numspecies * ones(size(communities, 2))
    @test subdiv(NormalisedAlpha(Metacommunity(communities)), [0])[!,
                                                                   :diversity] ≈
          numspecies * ones(size(communities, 2))
    qs = [0, 1, 2, Inf]
    sna = subdiv(NormalisedAlpha(Metacommunity(communities, Z1)), qs)
    @test nrow(sna) == length(qs) * size(communities, 2)
    for q in qs
        @test sna[isapprox.(sna[!, :q], q), :diversity] ≈
              ones(size(communities, 2))
    end
    @test subdiv(RawAlpha(Metacommunity(communities)), 0)[!, :diversity] ≈
          numspecies * vec(mapslices(v -> 1.0 / sum(v), communities; dims = 1))

    even = ones((numspecies, numcommunities)) / (numspecies * numcommunities)
    qs = [0, 1, 2, 3, 4, 5, 6, Inf]
    @test metadiv(NormalisedAlpha(Metacommunity(even)), qs)[!, :diversity] ≈
          numspecies * ones(length(qs))
    @test metadiv(RawAlpha(Metacommunity(even)), qs)[!, :diversity] ≈
          numspecies * numcommunities * ones(length(qs))
    md2 = metadiv(meta2, Inf)
    @test md2[md2[!, :measure] .== "Gamma", :diversity] ≈
          metadiv(Gamma(meta2), Inf)[!, :diversity]

    probs = reshape(mapslices(sum, communities; dims = 2),
                    size(communities, 1))
    @test metadiv(Gamma(Metacommunity(communities)), qs)[!, :diversity] ≈
          qD(probs, qs)
    @test metadiv(Gamma(Metacommunity(communities, Z1)), qs)[!, :diversity] ≈
          qDZ(probs, qs, Z1)

    Z = rand(numspecies, numspecies)
    @test metadiv(Gamma(Metacommunity(communities, Z)), qs)[!, :diversity] ≈
          qDZ(probs, qs, Z)

    colweights = rand(numcommunities)
    colweights /= sum(colweights)
    allthesame = probs * colweights'
    @test metadiv(RawBeta(Metacommunity(allthesame, Z)), qs)[!, :diversity] ≈
          1.0 ./ qD(colweights, 2 .- qs)
    @test metadiv(NormalisedBeta(Metacommunity(allthesame, Z)), qs)[!,
                                                                    :diversity] ≈
          ones(length(qs))
    @test metadiv(NormalisedRho(Metacommunity(allthesame, Z)), qs)[!,
                                                                   :diversity] ≈
          ones(length(qs))
    @test metadiv(RawRho(Metacommunity(allthesame, Z)), qs)[!, :diversity] ≈
          qD(colweights, qs)

    communitylist = rand(1:numcommunities, numspecies)
    distinct = zeros(Float64, (numspecies, numcommunities))
    for i in 1:numspecies
        distinct[i, communitylist[i]] = weights[i]
    end

    @test metadiv(RawRho(Metacommunity(distinct)), qs)[!, :diversity] ≈
          ones(length(qs))
    subnr = subdiv(NormalisedRho(Metacommunity(distinct)), qs)
    for q in qs
        @test subnr[isapprox.(subnr[!, :q], q), :diversity] ≈
              vec(sum(distinct; dims = 1))
    end
    @test metadiv(NormalisedBeta(Metacommunity(distinct)), qs)[!, :diversity] ≈
          qD(reshape(sum(distinct; dims = 1), numcommunities), qs)
    @test metadiv(RawBeta(Metacommunity(distinct)), qs)[!, :diversity] ≈
          ones(length(qs))

    # many (unexported!) diversity levels not yet implemented
    @test_throws ErrorException Diversity.communityDiversity(nab)
end

@testset "Plot recipe" begin
    # The recipe is defined on a *tuple*, so the measure and the order go in together.
    mc = Metacommunity(manyweights)
    @test plot((NormalisedRho(mc), 1)) isa Plots.Plot
    @test plot((Gamma(mc), 0)) isa Plots.Plot

    # And it needs coordinates. A partition with no spatial data makes some up, but only because
    # `getcoords` is wired to `coordinates` in `src/EcoBase.jl` — EcoBase's own fallback returns the
    # partition itself, which Plots cannot use.
    @test getcoords(getpartition(mc)) isa AbstractMatrix
    @test size(getcoords(getpartition(mc))) == (countsubcommunities(mc), 2)

    # The recipe plots *subcommunity* diversities against the partition's coordinates, so it needs
    # one value per subcommunity.
    @test nrow(subdiv(NormalisedRho(mc), 1)) == countsubcommunities(mc)
end

@testset "Choosing what comes back" begin
    # The result is built as columns and then materialised into whatever table the caller asks for.
    # A DataFrame stays the default, so every existing call is unaffected.
    mc = Metacommunity([0.1 0.2; 0.2 0.1; 0.2 0.2])
    a = NormalisedAlpha(mc)

    @test subdiv(a, 1) isa DataFrame
    @test subdiv(DataFrame, a, 1) == subdiv(a, 1)
    @test inddiv(DataFrame, a, 1) == inddiv(a, 1)
    @test metadiv(DataFrame, a, 1) == metadiv(a, 1)

    # Any Tables sink works, and agrees with the DataFrame column for column.
    ct = subdiv(Tables.columntable, a, [0, 1])
    df = subdiv(a, [0, 1])
    @test keys(ct) == Tuple(Symbol.(names(df)))
    @test all(collect(ct[k]) == df[!, k] for k in keys(ct))

    # What is handed to the sink is a Tables source in its own right, which is what lets
    # CSV.write and friends take a result with no conversion step.
    @test Tables.istable(typeof(subdiv(a, 1)))

    # The sink threads through the wrappers and through the combined entry point.
    @test norm_sub_alpha(DataFrame, mc, 1) == norm_sub_alpha(mc, 1)
    @test meta_gamma(DataFrame, mc, 1) == meta_gamma(mc, 1)
    levels = [subcommunityDiversity, metacommunityDiversity]
    @test diversity(DataFrame, levels, [ᾱ, Γ], mc, [0, 1]) ==
          diversity(levels, [ᾱ, Γ], mc, [0, 1])

    # The all-seven forms take a sink too, at every scale.
    for f in (inddiv, subdiv, metadiv)
        @test f(DataFrame, mc, 1) == f(mc, 1)
        @test collect(f(Tables.columntable, mc, [0, 1]).diversity) ==
              f(mc, [0, 1]).diversity
    end

    # Individual diversities are a level like any other, so `diversity` can ask for them.
    ind = diversity([individualDiversity], [ᾱ], mc, 1)
    @test nrow(ind) == counttypes(mc) * countsubcommunities(mc)
    @test ind == inddiv(ᾱ(mc), 1)
    @test collect(diversity(Tables.columntable, [individualDiversity], [ᾱ], mc,
                            1).diversity) == ind.diversity

    # Several measures, orders and levels in one call: each measure is built once and asked for
    # every level, which is the reason this entry point exists.
    combined = diversity(levels, [ᾱ, ρ̄, Γ], mc, [0, 1, 2])
    @test length(unique(combined.measure)) == 3
    @test length(unique(combined.q)) == 3
    @test length(unique(combined.partition_level)) == 2
end

@testset "Every wrapper takes a sink, and takes the right one" begin
    # The fourteen sink methods were generated mechanically, so a copy-paste error in one of them
    # would be invisible without checking all fourteen: each must reach the same measure and the
    # same scale as the sinkless method it shadows.
    species = ["a", "b", "c"]
    sites = ["north", "south"]
    mc = Metacommunity([0.1 0.2; 0.2 0.1; 0.2 0.2], UniqueTypes(species),
                       Subcommunities(sites))
    subs = [norm_sub_alpha, raw_sub_alpha, norm_sub_beta, raw_sub_beta,
        norm_sub_rho, raw_sub_rho, sub_gamma]
    metas = [norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta,
        norm_meta_rho, raw_meta_rho, meta_gamma]

    for f in vcat(subs, metas)
        @test f(DataFrame, mc, 1) == f(mc, 1)
        columns = f(Tables.columntable, mc, 1)
        @test collect(columns.diversity) == f(mc, 1).diversity
        @test collect(columns.measure) == f(mc, 1).measure
    end

    # A subcommunity wrapper must not have been wired to metadiv, or vice versa.
    @test all(f -> only(unique(f(mc, 1).partition_level)) == "subcommunity",
              subs)
    @test all(f -> only(unique(f(mc, 1).partition_level)) == "metacommunity",
              metas)
    # Seven measures, each named by one subcommunity and one metacommunity wrapper.
    @test length(unique(first(f(mc, 1).measure) for f in subs)) == 7
    @test [first(f(mc, 1).measure) for f in subs] ==
          [first(f(mc, 1).measure) for f in metas]
end

@testset "Result rows are in the documented order" begin
    # Individual diversities are one row per type per subcommunity, with the types cycling fastest
    # -- the order the old row-at-a-time construction produced, and what rdiversity expects.
    species = ["a", "b", "c"]
    sites = ["north", "south"]
    mc = Metacommunity([0.1 0.2; 0.2 0.1; 0.2 0.2], UniqueTypes(species),
                       Subcommunities(sites))
    ind = inddiv(ᾱ(mc), 1)
    @test nrow(ind) == length(species) * length(sites)
    @test ind.type_name == repeat(species, outer = length(sites))
    @test ind.partition_name == repeat(sites, inner = length(species))
    @test ind.diversity == vec(inddiv(ᾱ(mc), 1).diversity)

    # Subcommunity diversities are one row per subcommunity, in order, with no type named.
    sub = subdiv(ᾱ(mc), 1)
    @test sub.partition_name == sites
    @test all(isempty, sub.type_name)

    # Several orders are stacked in the order asked for, not interleaved.
    many = subdiv(ᾱ(mc), [0, 1, 2])
    @test many.q == repeat([0, 1, 2], inner = length(sites))
    @test many.partition_name == repeat(sites, outer = 3)

    # All seven measures come back in their long-standing order.
    @test unique(subdiv(mc, 1).measure) ==
          ["RawAlpha", "NormalisedAlpha", "RawBeta", "NormalisedBeta", "RawRho",
        "NormalisedRho", "Gamma"]

    # A level with no implementation is refused rather than silently skipped.
    @test_throws ErrorException diversity([Diversity.communityDiversity], [ᾱ],
                                          mc, 1)
end

end
