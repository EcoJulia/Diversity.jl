# SPDX-License-Identifier: BSD-2-Clause

module ValidateRCall_rdiversity
using Test
using Statistics

using Diversity
using Diversity.Ecology
using Diversity.ShortNames
using DataFrames
using Phylo
using PopGen
using RCall

# The output columns both packages promise. Neither package's column set contains the other's,
# which is why the subset assertions these replace could never pass, in either direction: `div_type`
# is this package's name for what rdiversity calls `dat_id` - and the vocabularies differ too, ours
# saying "Arbitrary Z" where rdiversity says "UserGenerated" - while rdiversity carries four columns
# recording how a similarity matrix was derived from distances (`normalised`, `k`, `max_d`,
# `transformation`) that this package does not. So the shared contract is what is worth asserting,
# and asserting it is not vacuous: it is exactly what lets a result move between the two packages.
const SHARED_COLUMNS = Set(["diversity", "measure", "partition_level",
                               "partition_name", "q", "type_level", "type_name"
                           ])

# Check a Diversity result and an rdiversity result against the columns they both promise.
# Column *names* only. The `measure` column's *values* differ between the packages for every
# measure - "NormalisedBeta" against "normalised beta", "Gamma" against "gamma" - which this does
# not test. That divergence is known and is a question for the two maintainers, not an oversight.
function _checkcolumns(jmd, rmd)
    jcols = Set(map(string, names(jmd)))
    rcols = Set(rcopy(rcall(:colnames, rmd)))
    @test SHARED_COLUMNS ⊆ jcols
    @test SHARED_COLUMNS ⊆ rcols
    @test setdiff(jcols, Set(["div_type"])) ⊆ rcols
    return nothing
end

# Create a temporary directory to work in
libdir = mktempdir()

# Environment variable to avoid boring R package builds
global mustCrossvalidate = haskey(ENV, "JULIA_MUST_CROSSVALIDATE") &&
                           ENV["JULIA_MUST_CROSSVALIDATE"] == "1"

global skipR = false
if !skipR && !rcopy(R"require(ape)")
    rcall(Symbol(".libPaths"), libdir)
    reval("install.packages(\"ape\", lib=\"$libdir\", " *
          "repos=\"http://cran.r-project.org\")")
    global skipR = !rcopy(R"require(ape, lib.loc=c(\"$libdir\", .libPaths()))") &&
                   !mustCrossvalidate
    skipR && @warn "ape R package not installed and would not install, " *
          "skipping R crossvalidation"
end

if !skipR && !rcopy(R"require(rdiversity)")
    rcall(Symbol(".libPaths"), libdir)
    reval("install.packages(\"rdiversity\", lib=\"$libdir\", " *
          "repos=\"http://cran.r-project.org\")")
    global skipR = !rcopy(R"require(rdiversity, lib.loc=c(\"$libdir\", .libPaths()))") &&
                   !mustCrossvalidate
    skipR &&
        @warn "rdiversity R package not installed and would not install, " *
              "skipping R crossvalidation"
end

if !skipR && !rcopy(R"require(vegan)")
    rcall(Symbol(".libPaths"), libdir)
    reval("install.packages(\"vegan\", lib=\"$libdir\", " *
          "repos=\"http://cran.r-project.org\")")
    global skipR = !rcopy(R"require(vegan, lib.loc=c(\"$libdir\", .libPaths()))") &&
                   !mustCrossvalidate
    skipR && @warn "vegan R package not installed and would not install, " *
          "skipping R crossvalidation"
end

# adiv (which imports ade4) is needed only for pairwisedissimilarity(), so
# failing to get it skips that comparison and not the rest
global skipadiv = skipR
if !skipadiv
    # adiv imports rgl, which must not try to open a display
    R"options(rgl.useNULL = TRUE)"
    if !rcopy(R"suppressMessages(requireNamespace(\"adiv\", quietly = TRUE))")
        rcall(Symbol(".libPaths"), libdir)
        reval("install.packages(\"adiv\", lib=\"$libdir\", " *
              "repos=\"http://cran.r-project.org\")")
        global skipadiv = !rcopy(R"suppressMessages(requireNamespace(\"adiv\", lib.loc=c(\"$libdir\", .libPaths()), quietly = TRUE))") &&
                          !mustCrossvalidate
        skipadiv && @warn "adiv R package not installed and would not " *
              "install, skipping its crossvalidation"
    end
end

if !skipR
    # vegan >= 2.7 changed method = "gower" to range-standardise columns via
    # decostand() and drop zero-range (tied / double-zero) columns from the
    # denominator, returning NA when all columns are tied. Diversity follows the
    # classic Gower (1971) convention, where tied columns contribute 0 but stay
    # in M (the total number of types). gowerOld() reproduces that for the two
    # subcommunities compared here: gower is the column-mean of the
    # range-standardised Manhattan distance, so each differing column
    # contributes 1 and each tied column (NaN -> 0) contributes 0.
    R"""
    gowerOld <- function(m) {
        xs <- vegan::decostand(m, "range", 2, na.rm = TRUE)
        xs[is.na(xs)] <- 0
        mean(abs(xs[1, ] - xs[2, ]))
    }
    """

    # Run diversity comparisons with vegan
    @testset "Two community diversity measures from vegan" begin
        @testset "Random pop $i" for i in 2 .^ (1:5)
            types = i
            sc = 2
            pops = rand(types, sc)
            pops ./= sum(pops)

            @rput pops
            R"""
            rp = t(pops)
            rj = vegdist(rp, method = "jaccard")[1]
            rj1 = vegdist(rp > 0, method = "jaccard")[1]
            rg = gowerOld(rp)
            rg1 = gowerOld(rp > 0)
            rag = vegdist(rp, method = "altGower")[1]
            rag1 = vegdist(rp > 0, method = "altGower")[1]
            """
            @rget rj
            @rget rj1
            @rget rg
            @rget rg1
            @rget rag
            @rget rag1

            jj = jaccard(pops).diversity[1]
            jj1 = jaccard([x > 0 ? 1 : 0 for x in pops]).diversity[1]
            jg = gower(pops, countzeros = true).diversity[1]
            jg1 = gower([x > 0 ? 1 : 0 for x in pops], countzeros = true).diversity[1]
            jag = gower(pops, countzeros = false).diversity[1]
            jag1 = gower([x > 0 ? 1 : 0 for x in pops], countzeros = false).diversity[1]

            @test jj ≈ 1.0 - rj
            @test jj1 ≈ 1.0 - rj1
            @test jg ≈ rg
            @test jg1 ≈ rg1
            @test jag ≈ rag
            @test jag1 ≈ rag1

            for j in 1:sc
                vals = pops[:, j] .< median(pops[:, j]) / 2
                for k in 1:types
                    if vals[k]
                        pops[k, j] = 0.0
                    end
                end
            end
            pops ./= sum(pops)

            @rput pops
            R"""
            rp = t(pops)
            rj = vegdist(rp, method = "jaccard")[1]
            rj1 = vegdist(rp > 0, method = "jaccard")[1]
            rg = gowerOld(rp)
            rg1 = gowerOld(rp > 0)
            rag = vegdist(rp, method = "altGower")[1]
            rag1 = vegdist(rp > 0, method = "altGower")[1]
            """
            @rget rj
            @rget rj1
            @rget rg
            @rget rg1
            @rget rag
            @rget rag1

            jj = jaccard(pops).diversity[1]
            jj1 = jaccard([x > 0 ? 1 : 0 for x in pops]).diversity[1]
            jg = gower(pops, countzeros = true).diversity[1]
            jg1 = gower([x > 0 ? 1 : 0 for x in pops], countzeros = true).diversity[1]
            jag = gower(pops, countzeros = false).diversity[1]
            jag1 = gower([x > 0 ? 1 : 0 for x in pops], countzeros = false).diversity[1]

            @test jj ≈ 1.0 - rj
            @test jj1 ≈ 1.0 - rj1
            @test jg ≈ rg
            @test jg1 ≈ rg1
            @test jag ≈ rag
            @test jag1 ≈ rag1
        end
    end

    # pairwisedissimilarity() against adiv's dsimcom() (the four forms of
    # Pavoine & Ricotta 2014), ade4's disc() (Rao's DISC), and vegan's vegdist()
    # for the classic indices they reduce to without a similarity between types
    if !skipadiv
        @testset "Pairwise dissimilarity from adiv, ade4 and vegan" begin
            @testset "Random pop $i" for i in 1:10
                types = rand(2:(i * 5))
                sc = rand(2:(i * 3 + 2))
                ab = rand(types, sc)
                # Make sure not to remove all of the non-zeros from any column
                for j in 1:sc
                    ab[ab[:, j] .< median(ab[:, j]) / 2, j] .= 0.0
                end
                ab ./= sum(ab)
                pa = Float64.(ab .> 0)
                pa ./= sum(pa)
                # dsimcom() needs a positive semi-definite similarity in
                # [0, 1]: the cosines between nonnegative feature vectors
                # A small feature of its own keeps every type distinct
                nfeatures = rand(1:(2 * types))
                F = hcat(0.1 .* ((1:types) .== (1:types)'),
                         rand(types, nfeatures) .*
                         (rand(types, nfeatures) .< 0.5))
                F ./= sqrt.(sum(abs2, F, dims = 2))
                Z = F * F'
                Z = (Z + Z') / 2
                for k in 1:types
                    Z[k, k] = 1.0
                end

                @rput ab pa Z
                R"""
                comm <- t(ab)
                dsim <- function(option) lapply(1:4, function(m)
                    unname(as.matrix(adiv::dsimcom(comm, Z, method = m,
                                                   option = option,
                                                   type = "dissimilarity"))))
                rabs <- dsim("absolute")
                rrel <- dsim("relative")
                rdisc <- unname(as.matrix(ade4::disc(as.data.frame(ab),
                                                     as.dist(sqrt(2 * (1 - Z))))))
                rbinj <- unname(as.matrix(vegdist(t(pa) > 0, method = "jaccard")))
                rbinb <- unname(as.matrix(vegdist(t(pa) > 0, method = "bray")))
                rhorn <- unname(as.matrix(vegdist(comm, method = "horn")))
                """
                @rget rabs rrel rdisc rbinj rbinb rhorn

                # dsimcom()'s methods 1 to 4
                forms = (:sokalsneath, :jaccard, :sorensen, :ochiai)
                dabs = pairwisedissimilarity(ab, Z)
                drel = pairwisedissimilarity(ab, Z; relative = true)
                for (m, form) in enumerate(forms)
                    @test dabs[form] ≈ rabs[m] atol = 1e-12
                    @test drel[form] ≈ rrel[m] atol = 1e-12
                end
                # disc() squares the distances it is given, so √(2(1 - Z))
                # makes its Δ twice ours, and it reports √DISC
                @test sqrt.(2 .* dabs.rao) ≈ rdisc atol = 1e-12
                # With no similarity between types: binary Jaccard and
                # Bray-Curtis (Sørensen), and Morisita-Horn (relative Sørensen)
                d01 = pairwisedissimilarity(pa)
                @test d01.jaccard ≈ rbinj atol = 1e-12
                @test d01.sorensen ≈ rbinb atol = 1e-12
                @test pairwisedissimilarity(ab; relative = true).sorensen ≈
                      rhorn atol = 1e-12
            end
        end
    end

    # Run diversity comparisons on increasing numbers of types and subcommunities
    @testset "RCall - testing boydorr/rdiversity" begin
        @testset "Random pop $i" for i in 1:10
            types = rand(2:(i * 10))
            sc = rand(2:(i * 10))
            pops = rand(types, sc)
            Zasym = rand(types, types)
            Zsym = Zasym / 2 .+ Zasym' / 2
            for j in 1:types
                Zsym[j, j] = 1.0
            end
            Zsym[Zsym .< median(Zsym)] .= 0.0
            # Make sure not to remove all of the non-zeros from any column
            for j in 1:sc
                pops[pops[:, j] .< median(pops[:, j]) / 2, j] .= 0.0
            end
            pops /= sum(pops)
            qs = sort([rand(7) * 10..., 0, 1, Inf])
            names = ["symmetric", "asymmetric"]
            Zs = [Zsym, Zasym]
            @testset "Z matrix - $(names[k])" for k in axes(Zs, 1)
                Z = Zs[k]
                meta = Metacommunity(pops, Z)
                diversities = Dict(:raw_alpha => α(meta),
                                   :norm_alpha => ᾱ(meta),
                                   :raw_beta => β(meta),
                                   :norm_beta => β̄(meta),
                                   :raw_rho => ρ(meta),
                                   :norm_rho => ρ̄(meta),
                                   :raw_gamma => Γ(meta))
                # Create the metacommunity in R
                sim = rcall(:similarity, Z)
                r_meta = rcall(:metacommunity, pops, sim)

                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta)
                    # Check the metacommunity diversity
                    @test metadiv(juliadiv, qs)[!, :diversity] ≈
                          rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                    # and subcommunity diversity
                    @test subdiv(juliadiv, qs)[!, :diversity] ≈
                          rcopy(rcall(:subdiv, r_div, qs)[:diversity])
                end
            end
        end

        @testset "Trying out empty types and subcommunities" begin
            @testset "Random pop $i" for i in 1:10
                types = 10
                sc = 10
                pops = rand(types, sc)
                Zasym = rand(types, types)
                Z = Zasym / 2 + Zasym' / 2
                for j in 1:types
                    Z[j, j] = 1.0
                end
                Z[Z .< median(Z)] .= 0.0
                # Make sure not to remove all of the non-zeros from any column
                for j in 1:sc
                    vals = pops[:, j] .< median(pops[:, j]) / 2
                    for k in 1:types
                        if vals[k]
                            pops[k, j] = 0.0
                        end
                    end
                end
                qs = sort([rand(7) * 10..., 0, 1, Inf])

                # Check they match when there's an empty type
                pops[rand(1:types), :] .= 0
                pops /= sum(pops)
                meta = Metacommunity(pops, Z)
                diversities = Dict(:raw_alpha => α(meta),
                                   :norm_alpha => ᾱ(meta),
                                   :raw_beta => β(meta),
                                   :norm_beta => β̄(meta),
                                   :raw_rho => ρ(meta),
                                   :norm_rho => ρ̄(meta),
                                   :raw_gamma => Γ(meta))
                # Create the metacommunity in R
                sim = rcall(:similarity, Z)
                r_meta = rcall(:metacommunity, pops, sim)

                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta)
                    # Check the metacommunity diversity
                    jmd = metadiv(juliadiv, qs)
                    rmd = rcall(:metadiv, r_div, qs)
                    _checkcolumns(jmd, rmd)

                    @test jmd[!, :diversity] ≈ rcopy(rmd[:diversity])
                    # and subcommunity diversity
                    @test subdiv(juliadiv, qs)[!, :diversity] ≈
                          rcopy(rcall(:subdiv, r_div, qs)[:diversity])
                end

                # Check they match when there's an empty subcommunity too
                pops[:, rand(1:sc)] .= 0
                pops /= sum(pops)
                meta = Metacommunity(pops, Z)
                diversities = Dict(:raw_alpha => α(meta),
                                   :norm_alpha => ᾱ(meta),
                                   :raw_beta => β(meta),
                                   :norm_beta => β̄(meta),
                                   :raw_rho => ρ(meta),
                                   :norm_rho => ρ̄(meta),
                                   :raw_gamma => Γ(meta))
                # Create the metacommunity in R
                sim = rcall(:similarity, Z)
                r_meta = rcall(:metacommunity, pops, sim)

                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta)
                    # Check out metacommunity diversity
                    @test metadiv(juliadiv, qs)[!, :diversity] ≈
                          rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                    # and subcommunity diversity
                    sdj = subdiv(juliadiv, qs)[!, :diversity]
                    sdr = rcopy(rcall(:subdiv, r_div, qs)[:diversity])
                    for (r, j) in zip(sdr, sdj)
                        @test isnan(j) == isnan(r) && (isnan(j) || j ≈ r)
                    end
                end
            end
        end
    end

    # Run phylogenetic comparisons
    @testset "RCall - testing Phylogenetics with boydorr/rdiversity" begin
        @testset "Random phylogeny $i" for i in 1:10
            types = rand(2:(i * 5))
            nu = Nonultrametric(types)
            tree = rand(nu)
            sc = rand(2:(i * 10))
            pops = rand(types, sc)
            for j in 1:sc
                pops[pops[:, j] .< median(pops[:, j]) / 2, j] .= 0.0
            end
            pops ./= sum(pops)
            qs = sort([rand(7) * 10..., 0, 1, Inf])
            meta = Metacommunity(pops, PhyloBranches(tree))
            diversities = Dict(:raw_alpha => α(meta),
                               :norm_alpha => ᾱ(meta),
                               :raw_beta => β(meta),
                               :norm_beta => β̄(meta),
                               :raw_rho => ρ(meta),
                               :norm_rho => ρ̄(meta),
                               :raw_gamma => Γ(meta))
            # Create the metacommunity in R
            @rput pops
            @rput tree
            tipnames = [getnodename(tree, node)
                        for node in traversal(tree, preorder)
                        if isleaf(tree, node)]
            @rput tipnames
            r_meta = R"""
            rownames(pops) <- tipnames
            metacommunity(pops, phy2branch(tree, pops))
            """

            for (r_func, juliadiv) in diversities
                r_div = rcall(r_func, r_meta)
                # Check the metacommunity diversity
                @test metadiv(juliadiv, qs)[!, :diversity] ≈
                      rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                # and subcommunity diversity
                @test subdiv(juliadiv, qs)[!, :diversity] ≈
                      rcopy(rcall(:subdiv, r_div, qs)[:diversity])
            end
        end

        @testset "Trying out empty types and subcommunities" begin
            types = 10
            nu = Nonultrametric(types)
            tree = rand(nu)
            sc = 10
            pops = rand(types, sc)
            qs = sort([rand(7) * 10..., 0, 1, Inf])

            # Check they match when there's an empty type
            pops[rand(1:types), :] .= 0
            pops /= sum(pops)
            meta = Metacommunity(pops, PhyloBranches(tree))
            diversities = Dict(:raw_alpha => α(meta),
                               :norm_alpha => ᾱ(meta),
                               :raw_beta => β(meta),
                               :norm_beta => β̄(meta),
                               :raw_rho => ρ(meta),
                               :norm_rho => ρ̄(meta),
                               :raw_gamma => Γ(meta))
            # Create the metacommunity in R
            @rput pops
            @rput tree
            tipnames = [getnodename(tree, node)
                        for node in traversal(tree, preorder)
                        if isleaf(tree, node)]
            @rput tipnames
            r_meta = R"""
            rownames(pops) <- tipnames
            r_meta <- metacommunity(pops, phy2branch(tree, pops))
            """

            for (r_func, juliadiv) in diversities
                r_div = rcall(r_func, r_meta)
                # Check the metacommunity diversity
                jmd = metadiv(juliadiv, qs)
                rmd = rcall(:metadiv, r_div, qs)
                _checkcolumns(jmd, rmd)

                @test jmd[!, :diversity] ≈ rcopy(rmd[:diversity])
                # and subcommunity diversity
                @test subdiv(juliadiv, qs)[!, :diversity] ≈
                      rcopy(rcall(:subdiv, r_div, qs)[:diversity])
            end

            # Check they match when there's an empty subcommunity too
            pops[:, rand(1:sc)] .= 0
            pops ./= sum(pops)
            meta = Metacommunity(pops, PhyloBranches(tree))
            diversities = Dict(:raw_alpha => α(meta),
                               :norm_alpha => ᾱ(meta),
                               :raw_beta => β(meta),
                               :norm_beta => β̄(meta),
                               :raw_rho => ρ(meta),
                               :norm_rho => ρ̄(meta),
                               :raw_gamma => Γ(meta))
            # Create the metacommunity in R
            @rput pops
            @rput tree
            tipnames = [getnodename(tree, node)
                        for node in traversal(tree, preorder)
                        if isleaf(tree, node)]
            @rput tipnames
            r_meta = R"""
            rownames(pops) <- tipnames
            metacommunity(pops, phy2branch(tree, pops))
            """
            for (r_func, juliadiv) in diversities
                r_div = rcall(r_func, r_meta)
                # Check out metacommunity diversity
                @test metadiv(juliadiv, qs)[!, :diversity] ≈
                      rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                # and subcommunity diversity
                sdj = subdiv(juliadiv, qs)[!, :diversity]
                sdr = rcopy(rcall(:subdiv, r_div, qs)[:diversity])
                for (r, j) in zip(sdr, sdj)
                    @test isnan(j) == isnan(r) && (isnan(j) || j ≈ r)
                end
            end
        end
    end

    # Run genetic (VCF) comparisons. For every VCF file in test/data, Julia reads
    # it with PopGen and builds a GeneticType with biallelic Manhattan distances,
    # which must match rdiversity's gen2dist(vcf, biallelic = TRUE) + dist2sim.
    # Drop additional .vcf files into test/data to cross-validate them too.
    @testset "RCall - testing Genetics with boydorr/rdiversity" begin
        datadir = joinpath(@__DIR__, "data")
        vcffiles = sort([f
                         for f in readdir(datadir; join = true)
                         if endswith(f, ".vcf")])
        @testset "VCF $(basename(file))" for file in vcffiles
            pd = PopGen.vcf(file; silent = true)
            gv = GeneticType(pd; distance = :manhattan)
            jnames = gettypenames(gv, true)

            # Build the matching R similarity from the *same* PopData that Julia
            # used, via vcf_dataframe (RCall converts the DataFrame to an R
            # data.frame that gen2dist accepts). It depends only on the genotypes.
            vcf = vcf_dataframe(pd)
            @rput vcf
            R"""
            gsim <- dist2sim(gen2dist(vcf, biallelic = TRUE), "linear")
            """

            # Cross-validate against a range of random abundances for this VCF
            @testset "Random abundances $j" for j in 1:5
                sc = rand(2:5)
                pops = rand(length(jnames), sc)
                pops ./= sum(pops)
                meta = Metacommunity(pops, gv)
                qs = sort([rand(7) * 10..., 0, 1, Inf])
                diversities = Dict(:raw_alpha => α(meta),
                                   :norm_alpha => ᾱ(meta),
                                   :raw_beta => β(meta),
                                   :norm_beta => β̄(meta),
                                   :raw_rho => ρ(meta),
                                   :norm_rho => ρ̄(meta),
                                   :raw_gamma => Γ(meta))

                @rput pops jnames
                r_meta = R"""
                rownames(pops) <- jnames
                metacommunity(pops, gsim)
                """

                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta)
                    # Check the metacommunity diversity
                    @test metadiv(juliadiv, qs)[!, :diversity] ≈
                          rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                    # and subcommunity diversity
                    @test subdiv(juliadiv, qs)[!, :diversity] ≈
                          rcopy(rcall(:subdiv, r_div, qs)[:diversity])
                end
            end
        end
    end
end
rm(libdir, force = true, recursive = true)

end
