module ValidateRCall_rdiversity
using Base.Test

using Diversity
using Diversity.ShortNames
using DataFrames

Rinstalled = false
try
    using RCall
    Rinstalled = true
catch
    warn("R not installed, skipping R cross-validation.")
end

if Rinstalled
    # Create a temporary directory to work in
    libdir = mktempdir();

    # Only run R on macs
    if is_apple()
        # Environment variable to avoid boring R package builds
        skipRinstall = haskey(ENV, "SKIP_R_INSTALL") && ENV["SKIP_R_INSTALL"] == "1"
        # Skip the (slow!) R package installation step
        if skipRinstall
            reval("library(rdiversity)");
        else
            rcall(Symbol(".libPaths"), libdir);
            rcall(Symbol("install.packages"),
                  ["devtools", "methods", "ggplot2", "ape",
                   "phangorn", "tidyr", "tibble", "phytools", "reshape2", "ggthemes"],
                  lib=libdir, repos="http://cran.r-project.org");
            reval("library(devtools, lib.loc=c(\"$libdir\", .libPaths()))");
            rcall(:install_git, "https://github.com/boydorr/rdiversity.git", lib=libdir);
            reval("library(rdiversity, lib.loc=c(\"$libdir\", .libPaths()))");
        end

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
                    diversities = Dict(:raw_alpha  => α(meta), :norm_alpha => ᾱ(meta),
                                       :raw_beta   => β(meta), :norm_beta  => β̄(meta),
                                       :raw_rho    => ρ(meta), :norm_rho   => ρ̄(meta),
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
                diversities = Dict(:raw_alpha  => α(meta), :norm_alpha => ᾱ(meta),
                                   :raw_beta   => β(meta), :norm_beta  => β̄(meta),
                                   :raw_rho    => ρ(meta), :norm_rho   => ρ̄(meta),
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
                
                # Check they match when there's an empty subcommunity too
                pops[:, rand(1:sc)] = 0
                pops /= sum(pops)
                meta = Metacommunity(pops, Z)
                diversities = Dict(:raw_alpha  => α(meta), :norm_alpha => ᾱ(meta),
                                   :raw_beta   => β(meta), :norm_beta  => β̄(meta),
                                   :raw_rho    => ρ(meta), :norm_rho   => ρ̄(meta),
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
    end
    rm(libdir, force=true, recursive=true);
end

end
