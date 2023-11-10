module ValidateRCall_rdiversity
using Test
using Statistics

using Diversity
using Diversity.Ecology
using Diversity.ShortNames
using DataFrames
using Phylo
using RCall

# Create a temporary directory to work in
libdir = mktempdir();

# Environment variable to avoid boring R package builds
global mustCrossvalidate = haskey(ENV, "JULIA_MUST_CROSSVALIDATE") && ENV["JULIA_MUST_CROSSVALIDATE"] == "1"

global skipR = false
if !skipR && !rcopy(R"require(ape)")
    rcall(Symbol(".libPaths"), libdir);
    reval("install.packages(\"ape\", lib=\"$libdir\", " *
          "repos=\"http://cran.r-project.org\")");
    global skipR =
        !rcopy(R"require(ape, lib.loc=c(\"$libdir\", .libPaths()))") &&
        !mustCrossvalidate;
    skipR && @warn "ape R package not installed and would not install, " *
        "skipping R crossvalidation"
end

if !skipR && !rcopy(R"require(rdiversity)")
    rcall(Symbol(".libPaths"), libdir);
    reval("install.packages(\"rdiversity\", lib=\"$libdir\", " *
          "repos=\"http://cran.r-project.org\")");
    global skipR =
        !rcopy(R"require(rdiversity, lib.loc=c(\"$libdir\", .libPaths()))") &&
        !mustCrossvalidate;
    skipR && @warn "rdiversity R package not installed and would not install, " *
        "skipping R crossvalidation"
end

if !skipR && !rcopy(R"require(vegan)")
    rcall(Symbol(".libPaths"), libdir);
    reval("install.packages(\"vegan\", lib=\"$libdir\", " *
          "repos=\"http://cran.r-project.org\")");
    global skipR =
        !rcopy(R"require(vegan, lib.loc=c(\"$libdir\", .libPaths()))") &&
        !mustCrossvalidate;
    skipR && @warn "vegan R package not installed and would not install, " *
        "skipping R crossvalidation"
end

if !skipR
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
            rg = vegdist(rp, method = "gower")[1]
            rg1 = vegdist(rp > 0, method = "gower")[1]
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
                vals = pops[:, j] .< median(pops[:, j])/2
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
            rg = vegdist(rp, method = "gower")[1]
            rg1 = vegdist(rp > 0, method = "gower")[1]
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

    # Run diversity comparisons on increasing numbers of types and subcommunities
    @testset "RCall - testing boydorr/rdiversity" begin
        @testset "Random pop $i" for i in 1:10
            types = rand(2:(i*10))
            sc = rand(2:(i*10))
            pops = rand(types, sc)
            Zasym = rand(types, types)
            Zsym = Zasym/2 .+ Zasym'/2
            for j in 1:types
                Zsym[j, j] = 1.0
            end
            Zsym[Zsym .< median(Zsym)] .= 0.0
            # Make sure not to remove all of the non-zeros from any column
            for j in 1:sc
                pops[pops[:, j] .< median(pops[:, j])/2, j] .= 0.0
            end
            pops /= sum(pops)
            qs = sort([rand(7)*10..., 0, 1, Inf])
            names = ["symmetric", "asymmetric"]
            Zs = [Zsym, Zasym]
            @testset "Z matrix - $(names[k])" for k in axes(Zs, 1)
                Z = Zs[k]
                meta = Metacommunity(pops, Z)
                diversities = Dict(:raw_alpha  => α(meta),
                                   :norm_alpha => ᾱ(meta),
                                   :raw_beta   => β(meta),
                                   :norm_beta  => β̄(meta),
                                   :raw_rho    => ρ(meta),
                                   :norm_rho   => ρ̄(meta),
                                   :raw_gamma  => Γ(meta));
                # Create the metacommunity in R
                sim = rcall(:similarity, Z)
                r_meta = rcall(:metacommunity, pops, sim)

                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta);
                    # Check the metacommunity diversity
                    @test metadiv(juliadiv, qs)[!,:diversity] ≈
                        rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                    # and subcommunity diversity
                    @test subdiv(juliadiv, qs)[!,:diversity] ≈
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
                Z = Zasym/2 + Zasym'/2
                for j in 1:types
                    Z[j, j] = 1.0
                end
                Z[Z .< median(Z)] .= 0.0
                # Make sure not to remove all of the non-zeros from any column
                for j in 1:sc
                    vals = pops[:, j] .< median(pops[:, j])/2
                    for k in 1:types
                        if vals[k]
                            pops[k, j] = 0.0
                        end
                    end
                end
                qs = sort([rand(7)*10..., 0, 1, Inf])

                # Check they match when there's an empty type
                pops[rand(1:types), :] .= 0
                pops /= sum(pops)
                meta = Metacommunity(pops, Z)
                diversities = Dict(:raw_alpha  => α(meta),
                                :norm_alpha => ᾱ(meta),
                                :raw_beta   => β(meta),
                                :norm_beta  => β̄(meta),
                                :raw_rho    => ρ(meta),
                                :norm_rho   => ρ̄(meta),
                                :raw_gamma  => Γ(meta));
                # Create the metacommunity in R
                sim = rcall(:similarity, Z)
                r_meta = rcall(:metacommunity, pops, sim)

                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta);
                    # Check the metacommunity diversity
                    jmd = metadiv(juliadiv, qs);
                    rmd = rcall(:metadiv, r_div, qs);
                    @test_skip Set(map(string, names(jmd))) ⊆
                        Set(rcopy(rcall(:colnames, rmd)))
                        

                    @test jmd[!,:diversity] ≈ rcopy(rmd[:diversity])
                    # and subcommunity diversity
                    @test subdiv(juliadiv, qs)[!,:diversity] ≈
                        rcopy(rcall(:subdiv, r_div, qs)[:diversity])
                end

                # Check they match when there's an empty subcommunity too
                pops[:, rand(1:sc)] .= 0
                pops /= sum(pops)
                meta = Metacommunity(pops, Z)
                diversities = Dict(:raw_alpha  => α(meta),
                                :norm_alpha => ᾱ(meta),
                                :raw_beta   => β(meta),
                                :norm_beta  => β̄(meta),
                                :raw_rho    => ρ(meta),
                                :norm_rho   => ρ̄(meta),
                                :raw_gamma  => Γ(meta));
                # Create the metacommunity in R
                sim = rcall(:similarity, Z)
                r_meta = rcall(:metacommunity, pops, sim)

                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta);
                    # Check out metacommunity diversity
                    @test metadiv(juliadiv, qs)[!,:diversity] ≈
                        rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                    # and subcommunity diversity
                    sdj = subdiv(juliadiv, qs)[!,:diversity]
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
            types = rand(2:(i*5))
            nu = Nonultrametric(types)
            tree = rand(nu)
            sc = rand(2:(i*10))
            pops = rand(types, sc)
            for j in 1:sc
                pops[pops[:, j] .< median(pops[:, j])/2, j] .= 0.0
            end
            pops ./= sum(pops)
            qs = sort([rand(7)*10..., 0, 1, Inf])
            meta = Metacommunity(pops, PhyloBranches(tree))
            diversities = Dict(:raw_alpha  => α(meta),
                               :norm_alpha => ᾱ(meta),
                               :raw_beta   => β(meta),
                               :norm_beta  => β̄(meta),
                               :raw_rho    => ρ(meta),
                               :norm_rho   => ρ̄(meta),
                               :raw_gamma  => Γ(meta));
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
                r_div = rcall(r_func, r_meta);
                # Check the metacommunity diversity
                @test metadiv(juliadiv, qs)[!,:diversity] ≈
                    rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                # and subcommunity diversity
                @test subdiv(juliadiv, qs)[!,:diversity] ≈
                    rcopy(rcall(:subdiv, r_div, qs)[:diversity])
            end
        end

        @testset "Trying out empty types and subcommunities" begin
            types = 10
            nu = Nonultrametric(types)
            tree = rand(nu)
            sc = 10
            pops = rand(types, sc)
            qs = sort([rand(7)*10..., 0, 1, Inf])

            # Check they match when there's an empty type
            pops[rand(1:types), :] .= 0
            pops /= sum(pops)
            meta = Metacommunity(pops, PhyloBranches(tree))
            diversities = Dict(:raw_alpha  => α(meta),
                               :norm_alpha => ᾱ(meta),
                               :raw_beta   => β(meta),
                               :norm_beta  => β̄(meta),
                               :raw_rho    => ρ(meta),
                               :norm_rho   => ρ̄(meta),
                               :raw_gamma  => Γ(meta));
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
                r_div = rcall(r_func, r_meta);
                # Check the metacommunity diversity
                jmd = metadiv(juliadiv, qs);
                rmd = rcall(:metadiv, r_div, qs);
                @test_skip Set(rcopy(rcall(:colnames, rmd))) ⊆
                    Set(map(string, names(jmd)))

                @test jmd[!,:diversity] ≈ rcopy(rmd[:diversity])
                # and subcommunity diversity
                @test subdiv(juliadiv, qs)[!,:diversity] ≈
                    rcopy(rcall(:subdiv, r_div, qs)[:diversity])
            end

            # Check they match when there's an empty subcommunity too
            pops[:, rand(1:sc)] .= 0
            pops ./= sum(pops)
            meta = Metacommunity(pops, PhyloBranches(tree))
            diversities = Dict(:raw_alpha  => α(meta),
                               :norm_alpha => ᾱ(meta),
                               :raw_beta   => β(meta),
                               :norm_beta  => β̄(meta),
                               :raw_rho    => ρ(meta),
                               :norm_rho   => ρ̄(meta),
                               :raw_gamma  => Γ(meta));
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
                r_div = rcall(r_func, r_meta);
                # Check out metacommunity diversity
                @test metadiv(juliadiv, qs)[!,:diversity] ≈
                    rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                # and subcommunity diversity
                sdj = subdiv(juliadiv, qs)[!,:diversity]
                sdr = rcopy(rcall(:subdiv, r_div, qs)[:diversity])
                for (r, j) in zip(sdr, sdj)
                    @test isnan(j) == isnan(r) && (isnan(j) || j ≈ r)
                end
            end
        end
    end
end
rm(libdir, force=true, recursive=true);

end
