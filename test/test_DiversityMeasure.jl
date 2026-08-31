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
    # These are the paper's own descriptions of the measures, at subcommunity level - which is
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

    # The descriptive aliases are the names the papers use, and are the same types - so a measure
    # reached through one spelling must be identical to the same measure reached through another.
    @test Diversity.Distinctiveness ≡ RawBeta
    @test Diversity.Redundancy ≡ RawRho
    @test Diversity.Representativeness ≡ NormalisedRho
    @test getFullName(Diversity.Representativeness(meta)) ==
          "representativeness"

    # `getASCIIName` strips the module prefix and the type parameters, so it names the *measure*
    # rather than the concrete parameterisation - which is what the output DataFrame carries.
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

@testset "Individual diversities are held as a rule, not an array" begin
    # The seven measures no longer materialise their ntypes x nplaces individual diversities;
    # each holds a closure over the arrays it needs. The values must be exactly what the
    # broadcast expressions they replaced produced, so those are written out here in full rather
    # than derived from the package -- a formula copied out of the source would only prove that
    # the source agrees with itself.
    for mc in (meta2, Metacommunity(ab3, GeneralTypes(Matrix(1.0I, 2, 2)), sc),
        Metacommunity(pop, ms, oc))
        zp = getordinariness!(mc)
        zP = getmetaordinariness!(mc)
        w = getweight(mc)
        ab = getabundance(mc)
        expected = Dict(RawAlpha => zp .^ -1,
                        NormalisedAlpha => w' ./ zp,
                        RawBeta => zp ./ zP,
                        NormalisedBeta => zp ./ (zP .* w'),
                        RawRho => zP ./ zp,
                        NormalisedRho => (zP .* w') ./ zp,
                        Gamma => fill!(similar(w), 1)' ./ zP)
        for (measure, want) in expected
            raw = Diversity.inddiv_raw(measure(mc), 1)
            @test raw isa Diversity.IndividualDiversities
            @test raw isa AbstractMatrix{Float64}
            @test size(raw) == size(ab)
            @test eltype(raw) == eltype(ab)
            @test Base.IndexStyle(typeof(raw)) == IndexCartesian()
            # Elementwise, and as a whole: `collect` and broadcasting are the two ways anything
            # downstream touches it, and both go through `getindex`.
            @test collect(raw) == want
            @test raw .* 2 == want .* 2
            @test all(raw[i, j] == want[i, j]
                      for i in axes(want, 1), j in axes(want, 2))
        end
    end

    # The whole point: building a measure over a metacommunity whose ordinariness is already
    # cached allocates nothing of that size. The weights vector is the only array involved.
    big = Metacommunity(fill(1 / 2000, 20, 100))
    getordinariness!(big)
    getmetaordinariness!(big)
    for measure in (RawAlpha, NormalisedAlpha, RawBeta, NormalisedBeta,
        RawRho, NormalisedRho, Gamma)
        measure(big)                       # compile it before counting bytes
        @test @allocated(measure(big)) < 8 * 20 * 100
    end
end

# A metacommunity that recomputes its ordinariness on every call and counts them. `SubAssemblage`,
# which `view` returns, is exactly this shape -- it inherits the generic `_getordinariness!` and so
# caches nothing -- but a counter is needed to assert how often it is asked.
const ORDCALLS = Ref(0)

mutable struct Uncached{FP, A, T, P} <:
               Diversity.API.AbstractMetacommunity{FP, A, A, T, P}
    types::T
    part::P
    ab::A
end

Diversity.API._gettypes(mc::Uncached) = mc.types
Diversity.API._getpartition(mc::Uncached) = mc.part
Diversity.API._getabundance(mc::Uncached, ::Bool) = mc.ab
function Diversity.API._getordinariness!(mc::Uncached)
    ORDCALLS[] += 1
    return Diversity.API._calcordinariness(mc.types, mc.ab, 1)
end

@testset "A measure asks for the ordinariness once, at construction" begin
    # The individual diversities being a closure over the ordinariness makes this worth pinning: a
    # rule that reached back into the metacommunity on every element would be catastrophic for a
    # subtype that does not cache, and nothing else in the suite would notice, since the package's
    # own `Metacommunity` does cache.
    Z = [1.0 0.5 0.0; 0.5 1.0 0.5; 0.0 0.5 1.0]
    p = [0.1 0.2; 0.2 0.1; 0.2 0.2]
    types, part = GeneralTypes(Z), Subcommunities(2)
    uncached = Uncached{Float64, typeof(p), typeof(types), typeof(part)}(types,
                                                                         part,
                                                                         p)
    cached = Metacommunity(p, GeneralTypes(Z), Subcommunities(2))

    for measure in (RawAlpha, NormalisedAlpha, RawBeta, NormalisedBeta,
        RawRho, NormalisedRho, Gamma)
        ORDCALLS[] = 0
        dm = measure(uncached)
        # Alpha and gamma need one of the two ordinarinesses, the betas and rhos both -- and
        # `_getmetaordinariness!` is itself defined in terms of `_getordinariness!`, hence two.
        built = ORDCALLS[]
        @test built in (1, 2)

        sub = subdiv(dm, [0, 1, 2])
        meta = metadiv(dm, [0, 1, 2])
        ind = inddiv(dm, 1)
        # Six passes over the individual diversities, and not one of them goes back to the
        # metacommunity.
        @test ORDCALLS[] == built

        reference = measure(cached)
        @test sub.diversity ≈ subdiv(reference, [0, 1, 2]).diversity
        @test meta.diversity ≈ metadiv(reference, [0, 1, 2]).diversity
        @test ind.diversity ≈ inddiv(reference, 1).diversity
    end

    # And the in-repo case that has no cache at all.
    whole = view(cached, sites = [1, 2])
    @test subdiv(NormalisedAlpha(whole), 1).diversity ≈
          subdiv(NormalisedAlpha(cached), 1).diversity
end

@testset "Repeated result columns are rules, not arrays" begin
    # Seven of the eight columns are a constant or a cycled short list. They are now held as rules,
    # so what has to be pinned is that they still produce exactly the `fill` and `repeat`
    # expressions they replaced -- written out here rather than derived from the package.
    lazymc = Metacommunity([0.1 0.2; 0.2 0.1; 0.2 0.2],
                           UniqueTypes(["a", "b", "c"]),
                           Subcommunities(["x", "y"]))
    lazydm = NormalisedAlpha(lazymc)
    nt, ns = 3, 2
    n = nt * ns
    cols = Diversity._inddiv_columns(lazydm, 1)

    @test collect(cols.div_type) == fill(getdiversityname(lazydm), n)
    @test collect(cols.measure) == fill("NormalisedAlpha", n)
    @test collect(cols.q) == fill(1, n)
    @test collect(cols.type_level) == fill("type", n)
    @test collect(cols.type_name) == repeat(["a", "b", "c"], outer = ns)
    @test collect(cols.partition_level) == fill("subcommunity", n)
    @test collect(cols.partition_name) == repeat(["x", "y"], inner = nt)
    @test cols.diversity isa Vector{Float64}

    @test cols.div_type isa Diversity.ConstantColumn
    @test cols.type_name isa Diversity.RepeatedColumn
    for col in (cols.div_type, cols.type_name, cols.partition_name)
        @test col isa AbstractVector
        @test length(col) == n
        @test size(col) == (n,)
        @test Base.IndexStyle(typeof(col)) == IndexLinear()
    end
    @test eltype(cols.type_name) == String
    @test eltype(cols.q) == Int

    # A subcommunity result has no type name to cycle, so all six of its repeated columns are
    # constant; only the subcommunity names and the diversities are real.
    subcols = Diversity._subdiv_columns(lazydm, 1)
    @test subcols.type_name isa Diversity.ConstantColumn
    @test collect(subcols.type_name) == ["", ""]
    @test subcols.partition_name == ["x", "y"]

    # A rule must not alias the names it was built from -- that is what the copy in the inner
    # constructor is for, and it is invisible without a test.
    names = ["a", "b", "c"]
    col = Diversity.RepeatedColumn(names, 1, 6)
    names[1] = "changed"
    @test col[1] == "a"
    @test collect(col) == repeat(["a", "b", "c"], outer = 2)
end

@testset "The DataFrame a caller gets back is unchanged" begin
    # The whole design rests on the sink deciding: DataFrame copies its columns by default, so the
    # rules are materialised back into ordinary mutable Vectors and a caller cannot tell. If that
    # ever stops being true, this is where it shows.
    lazymc = Metacommunity([0.1 0.2; 0.2 0.1; 0.2 0.2],
                           UniqueTypes(["a", "b", "c"]),
                           Subcommunities(["x", "y"]))
    for df in (inddiv(NormalisedAlpha(lazymc), 1),
        subdiv(NormalisedAlpha(lazymc), [0, 1]),
        metadiv(NormalisedAlpha(lazymc), 1))
        @test all(col -> col isa Vector, eachcol(df))
        @test eltype(df.measure) == String
        @test eltype(df.diversity) == Float64
    end

    df = inddiv(NormalisedAlpha(lazymc), 1)
    @test df.measure isa Vector{String}
    df.measure[1] = "mutated"
    @test df.measure[1] == "mutated"
    @test df.measure[2] == "NormalisedAlpha"    # and only the one row moved
    df.type_name[2] = "renamed"
    @test df.type_name[1] == "a"
    @test df.type_name[2] == "renamed"
end

@testset "Several orders or levels chain rather than concatenate" begin
    # Asking for more than one order, measure or level builds one part per combination. Joining them
    # with `vcat` materialised every rule the parts held, exactly when the result is largest, so
    # they are chained instead. What has to be pinned is that chaining is indistinguishable from the
    # concatenation it replaced -- including at the part boundaries, which is where an off-by-one
    # would hide.
    chainmc = Metacommunity([0.1 0.2; 0.2 0.1; 0.2 0.2],
                            UniqueTypes(["a", "b", "c"]),
                            Subcommunities(["x", "y"]))
    chaindm = NormalisedAlpha(chainmc)

    parts = [Diversity._subdiv_columns(chaindm, q) for q in [0, 1, 2]]
    for k in keys(first(parts))
        chained = Diversity._chaincolumn([part[k] for part in parts])
        @test chained == reduce(vcat, (part[k] for part in parts))
        @test length(chained) == sum(length(part[k]) for part in parts)
        @test eltype(chained) == eltype(first(parts)[k])
        @test chained isa Diversity.ChainedColumn
        @test Base.IndexStyle(typeof(chained)) == IndexLinear()
        # Every row individually, so a boundary cannot be papered over by a whole-vector compare.
        want = reduce(vcat, (part[k] for part in parts))
        @test all(chained[i] == want[i] for i in eachindex(want))
    end

    # Parts of different lengths -- a subcommunity result has one row per subcommunity, a
    # metacommunity result exactly one -- so the cumulative bounds have to be right, not assumed
    # uniform.
    mixedlengths = [Diversity._subdiv_columns(chaindm, 1),
        Diversity._metadiv_columns(chaindm, 1)]
    joined = Diversity._chaincolumn([part[:partition_name]
                                     for part in mixedlengths])
    @test joined == ["x", "y", ""]
    @test length(joined) == 3

    # Where the parts do not share a concrete type there is nothing to gain, so it falls back to
    # copying. Individual results cycle their type names where subcommunity results repeat one, so
    # asking for both is the case that reaches it -- and it must still be correct.
    mixedtypes = [Diversity._inddiv_columns(chaindm, 1),
        Diversity._subdiv_columns(chaindm, 1)]
    fellback = Diversity._chaincolumn([part[:type_name]
                                       for part in mixedtypes])
    @test !(fellback isa Diversity.ChainedColumn)
    @test fellback == vcat(repeat(["a", "b", "c"], outer = 2), ["", ""])

    # And end to end: the whole result, over two levels and three orders, is what it was. Note the
    # ordering `diversity` produces -- level outer, order inner, so every subcommunity row for every
    # order comes before the first metacommunity row.
    levels = [subcommunityDiversity, metacommunityDiversity]
    combined = diversity(levels, [NormalisedAlpha], chainmc, [0, 1, 2])
    separate = vcat(subdiv(NormalisedAlpha(chainmc), [0, 1, 2]),
                    metadiv(NormalisedAlpha(chainmc), [0, 1, 2]))
    @test size(combined) == size(separate)
    @test combined.diversity ≈ separate.diversity
    @test combined.partition_name == separate.partition_name
    @test combined.q == separate.q
    @test all(col -> col isa Vector, eachcol(combined))
end

@testset "Empty subcommunities change nothing" begin
    # A subcommunity with no individuals -- a sea cell in a species grid, an inactive cell in a
    # landscape simulation -- must not move any other subcommunity's answer, nor the
    # metacommunity's. This is what makes it safe to recognise them from their cached weight and
    # skip them, and it follows from the framework rather than the implementation: abundances are
    # relative to the whole metacommunity so zeros do not move the total, a subcommunity's measures
    # depend only on its own composition and the metacommunity as a whole, and a power mean ignores
    # zero-weight entries.
    occupied = rand(4, 3) .+ 0.1
    padded = zeros(4, 7)
    padded[:, [2, 4, 6]] .= occupied
    padded ./= sum(padded)
    occupied = padded[:, [2, 4, 6]]
    dead = [1, 3, 5, 7]
    Zmat = Matrix(1.0I, 4, 4)

    for (whole, part) in ((Metacommunity(padded), Metacommunity(occupied)),
                          (Metacommunity(padded, GeneralTypes(Zmat)),
                           Metacommunity(occupied, GeneralTypes(Zmat))))
        for measure in (RawAlpha, NormalisedAlpha, RawBeta, NormalisedBeta,
            RawRho, NormalisedRho, Gamma),
            q in (0, 0.5, 1, 2, Inf)
            full = subdiv(measure(whole), q).diversity
            @test full[[2, 4, 6]] ≈ subdiv(measure(part), q).diversity
            @test all(isnan, full[dead])
            @test metadiv(measure(whole), q).diversity[1] ≈
                  metadiv(measure(part), q).diversity[1]
        end
    end
end

@testset "Plot recipe" begin
    # The recipe is defined on a *tuple*, so the measure and the order go in together.
    mc = Metacommunity(manyweights)
    @test plot((NormalisedRho(mc), 1)) isa Plots.Plot
    @test plot((Gamma(mc), 0)) isa Plots.Plot

    # And it needs coordinates. A partition with no spatial data makes some up, but only because
    # `getcoords` is wired to `coordinates` in `src/EcoBase.jl` - EcoBase's own fallback returns the
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
    @test diversity(DataFrame, levels, [ᾱ, Γ], mc, [0, 1]) ==
          diversity(levels, [ᾱ, Γ], mc, [0, 1])

    # The all-seven forms take a sink too, at every scale.
    for f in (inddiv, subdiv, metadiv)
        @test f(DataFrame, mc, 1) == f(mc, 1)
        @test collect(f(Tables.columntable, mc, [0, 1]).diversity) ==
              f(mc, [0, 1]).diversity
    end

    # Individual diversities are a level like any other, so `diversity` can ask for them.
    ind = diversity([individualDiversity], [ᾱ], mc, 1)
    @test nrow(ind) == counttypes(mc) * countsubcommunities(mc)
    @test ind == inddiv(ᾱ(mc), 1)
    @test collect(diversity(Tables.columntable, [individualDiversity], [ᾱ], mc,
                            1).diversity) == ind.diversity

    # Several measures, orders and levels in one call: each measure is built once and asked for
    # every level, which is the reason this entry point exists.
    combined = diversity(levels, [ᾱ, ρ̄, Γ], mc, [0, 1, 2])
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
    ind = inddiv(ᾱ(mc), 1)
    @test nrow(ind) == length(species) * length(sites)
    @test ind.type_name == repeat(species, outer = length(sites))
    @test ind.partition_name == repeat(sites, inner = length(species))
    @test ind.diversity == vec(inddiv(ᾱ(mc), 1).diversity)

    # Subcommunity diversities are one row per subcommunity, in order, with no type named.
    sub = subdiv(ᾱ(mc), 1)
    @test sub.partition_name == sites
    @test all(isempty, sub.type_name)

    # Several orders are stacked in the order asked for, not interleaved.
    many = subdiv(ᾱ(mc), [0, 1, 2])
    @test many.q == repeat([0, 1, 2], inner = length(sites))
    @test many.partition_name == repeat(sites, outer = 3)

    # All seven measures come back in their long-standing order.
    @test unique(subdiv(mc, 1).measure) ==
          ["RawAlpha", "NormalisedAlpha", "RawBeta", "NormalisedBeta", "RawRho",
        "NormalisedRho", "Gamma"]

    # A level with no implementation is refused rather than silently skipped.
    @test_throws ErrorException diversity([Diversity.communityDiversity], [ᾱ],
                                          mc, 1)
end

end
