module ValidateRCall_rcall
macro R_str(script)
    script
end
macro rput(script)
    script
end
macro rget(script)
    script
end
export @R_str, @rput, @rget
end

module ValidateRCall_rdiversity
using Compat.Test
using Compat: @warn

using Diversity
using Diversity.ShortNames
using Diversity.Phylogenetics
using DataFrames
using Phylo

# Environment variable to avoid boring R package builds
mustCrossvalidate = haskey(ENV, "JULIA_MUST_CROSSVALIDATE") && ENV["JULIA_MUST_CROSSVALIDATE"] == "1"

# Only run R on unix or when R is installed because JULIA_MUST_CROSSVALIDATE is set to 1
skipR = !mustCrossvalidate && !is_unix()
Rinstalled = false
try
    skipR && error("Skipping R testing...")
    using RCall
    # Pkg.dir("Phylo") doesn't guarantee to give the correct directory for the package,
    # but nor does anything else insanely! Skips cross-validation on failure.
    include(joinpath(dirname(dirname(dirname(@__FILE__))), "Phylo", "src", "rcall.jl"))
    Rinstalled = true
catch
    if mustCrossvalidate
        error("R not installed, but JULIA_MUST_CROSSVALIDATE is set")
    else
        @warn "R or appropriate Phylo package not installed, skipping R cross-validation."
    end
    using ValidateRCall_rcall
end

if Rinstalled
    # Create a temporary directory to work in
    libdir = mktempdir();

    if !skipR && !rcopy(R"require(ape)")
        rcall(Symbol(".libPaths"), libdir);
        reval("install.packages(\"ape\", lib=\"$libdir\", " *
              "repos=\"http://cran.r-project.org\")");
        skipR = !rcopy(R"require(ape, lib.loc=c(\"$libdir\", .libPaths()))") &&
            !mustCrossvalidate;
        skipR && @warn "ape R package not installed and would not install, " *
            "skipping R crossvalidation"
    end

    if !skipR && !rcopy(R"require(rdiversity)")
        rcall(Symbol(".libPaths"), libdir);
        reval("install.packages(\"rdiversity\", lib=\"$libdir\", " *
              "repos=\"http://cran.r-project.org\")");
        skipR = !rcopy(R"require(rdiversity, lib.loc=c(\"$libdir\", .libPaths()))") &&
            !mustCrossvalidate;
        skipR && @warn "rdiversity R package not installed and would not install, " *
            "skipping R crossvalidation"
    end

    if !skipR
        # Run diversity comparisons on increasing numbers of types and subcommunities
        @testset "RCall - testing boydorr/rdiversity" begin
            @testset "Random rdiversity $i" for i in 1:20
                types = rand(2:(i*10))
                sc = rand(2:(i*10))
                pops = rand(types, sc)
                Zasym = rand(types, types)
                Zsym = Zasym/2 + Zasym'/2
                for j in 1:types
                    Zsym[j, j] = 1.0
                end
                Zsym[Zsym .< median(Zsym)] = 0.0
                # Make sure not to remove all of the non-zeros from any column
                for j in 1:sc
                    pops[pops[:, j] .< median(pops[:, j])/2, j] = 0.0
                end
                pops /= sum(pops)
                qs = sort([rand(7)*10..., 0, 1, Inf])
                names = ["symmetric", "asymmetric"]
                Zs = [Zsym, Zasym]
                @testset "Z matrix - $(names[k])" for k in 1:length(Zs)
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
                    r_meta = rcall(:metacommunity, pops, Z)

                    for (r_func, juliadiv) in diversities
                        r_div = rcall(r_func, r_meta);
                        # Check the metacommunity diversity
                        @test metadiv(juliadiv, qs)[:diversity] ≈
                            rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                        # and subcommunity diversity
                        @test subdiv(juliadiv, qs)[:diversity] ≈
                            rcopy(rcall(:subdiv, r_div, qs)[:diversity])
                    end
                end
            end

            @testset "Trying out empty types and subcommunities" begin
                types = 10
                sc = 10
                pops = rand(types, sc)
                Zasym = rand(types, types)
                Z = Zasym/2 + Zasym'/2
                for j in 1:types
                    Z[j, j] = 1.0
                end
                Z[Z .< median(Z)] = 0.0
                # Make sure not to remove all of the non-zeros from any column
                for j in 1:sc
                    pops[pops[:, j] .< median(pops[:, j])/2, j] = 0.0
                end
                qs = sort([rand(7)*10..., 0, 1, Inf])

                # Check they match when there's an empty type
                pops[rand(1:types), :] = 0
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
                r_meta = rcall(:metacommunity, pops, Z)

                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta);
                    # Check the metacommunity diversity
                    jmd = metadiv(juliadiv, qs);
                    rmd = rcall(:metadiv, r_div, qs);
                    @test Set(map(string, names(jmd))) ==
                        Set(rcopy(rcall(:colnames, rmd)))
                    @test jmd[:diversity] ≈ rcopy(rmd[:diversity])
                    # and subcommunity diversity
                    @test subdiv(juliadiv, qs)[:diversity] ≈
                        rcopy(rcall(:subdiv, r_div, qs)[:diversity])
                end

                # Check they match when there's an empty subcommunity too
                pops[:, rand(1:sc)] = 0
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
                r_meta = rcall(:metacommunity, pops, Z)

                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta);
                    # Check out metacommunity diversity
                    @test metadiv(juliadiv, qs)[:diversity] ≈
                        rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                    # and subcommunity diversity
                    sdj = subdiv(juliadiv, qs)[:diversity]
                    sdr = rcopy(rcall(:subdiv, r_div, qs)[:diversity])
                    for (r, j) in zip(sdr, sdj)
                        @test isnan(j) == isnan(r) && (isnan(j) || j ≈ r)
                    end
                end
            end
        end

        # Run phylogenetic comparisons
        @testset "RCall - testing .Phylogenetics with boydorr/rdiversity" begin
            @testset "Random phylogeny $i" for i in 1:10
                types = rand(2:(i*5))
                nu = Nonultrametric(types)
                tree = rand(nu)
                sc = rand(2:(i*10))
                pops = rand(types, sc)
                for j in 1:sc
                    pops[pops[:, j] .< median(pops[:, j])/2, j] = 0.0
                end
                pops /= sum(pops)
                qs = sort([rand(7)*10..., 0, 1, Inf])
                meta = Metacommunity(pops, PhyloTypes(tree))
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
                tipnames = getleafnames(tree)
                @rput tipnames
                r_meta = R"""
                rownames(pops) <- tipnames
                metacommunity(pops, tree)
                """

                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta);
                    # Check the metacommunity diversity
                    @test metadiv(juliadiv, qs)[:diversity] ≈
                        rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                    # and subcommunity diversity
                    @test subdiv(juliadiv, qs)[:diversity] ≈
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
                pops[rand(1:types), :] = 0
                pops /= sum(pops)
                meta = Metacommunity(pops, PhyloTypes(tree))
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
                tipnames = getleafnames(tree)
                @rput tipnames
                r_meta = R"""
                rownames(pops) <- tipnames
                r_meta <- metacommunity(pops, tree)
                """

                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta);
                    # Check the metacommunity diversity
                    jmd = metadiv(juliadiv, qs);
                    rmd = rcall(:metadiv, r_div, qs);
                    @test Set(map(string, names(jmd))) ==
                        Set(rcopy(rcall(:colnames, rmd)))
                    @test jmd[:diversity] ≈ rcopy(rmd[:diversity])
                    # and subcommunity diversity
                    @test subdiv(juliadiv, qs)[:diversity] ≈
                        rcopy(rcall(:subdiv, r_div, qs)[:diversity])
                end

                # Check they match when there's an empty subcommunity too
                pops[:, rand(1:sc)] = 0
                pops /= sum(pops)
                meta = Metacommunity(pops, PhyloTypes(tree))
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
                tipnames = getleafnames(tree)
                @rput tipnames
                r_meta = R"""
                rownames(pops) <- tipnames
                metacommunity(pops, tree)
                """
                for (r_func, juliadiv) in diversities
                    r_div = rcall(r_func, r_meta);
                    # Check out metacommunity diversity
                    @test metadiv(juliadiv, qs)[:diversity] ≈
                        rcopy(rcall(:metadiv, r_div, qs)[:diversity])
                    # and subcommunity diversity
                    sdj = subdiv(juliadiv, qs)[:diversity]
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

end
