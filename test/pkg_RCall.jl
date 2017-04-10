module ValidateRCall_rdiversity
using Base.Test

using Diversity
using Diversity.ShortNames
using RCall
using DataFrames

# We can't do the R testing on windows or macs because I can't install the packages...
is_me = haskey(ENV, "USER") && ENV["USER"] == "richardr"
if is_linux() || is_me
    if is_me
        R"library(rdiversity)";
    else
        libdir = mktempdir();
        rcall(Symbol(".libPaths"), libdir);
        rcall(Symbol("install.packages"),
              ["devtools", "methods", "ggplot2", "ape",
               "phangorn", "tidyr", "tibble", "phytools", "reshape2"],
              repos="http://cran.r-project.org", lib=libdir);
        @rput libdir
        R"library(devtools, lib.loc=c(libdir, .libPaths()))";
        rcall(Symbol("install_git"), "https://github.com/boydorr/rdiversity.git",
              lib=libdir);
        R"library(rdiversity, lib.loc=c(libdir, .libPaths()))";
    end

@testset "RCall - rdiversity" begin
    @testset "Random rdiversity $i" for i in 1:20
        types = rand(1:(i*10))
        sc = rand(1:(i*10))
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
        @testset "Z matrix ($names[k])" for k in 1:length(Zs)
            Z = Zs[k]
            meta = Metacommunity(pops, Z)
            ra = α(meta)
            na = ᾱ(meta)
            rb = β(meta)
            nb = β̄(meta)
            rr = ρ(meta)
            nr = ρ̄(meta)
            g  = Γ(meta)
            @rput pops
            @rput qs
            @rput Z
            # Checking out metacommunity diversities
            R"r_meta = metacommunity(pops, Z)
              r_rma = raw_meta_alpha(r_meta, qs)
              r_nma = norm_meta_alpha(r_meta, qs)
              r_rmb = raw_meta_beta(r_meta, qs)
              r_nmb = norm_meta_beta(r_meta, qs)
              r_rmr = raw_meta_rho(r_meta, qs)
              r_nmr = norm_meta_rho(r_meta, qs)
              r_mg  = meta_gamma(r_meta, qs)"
            @rget r_rma
            @test metadiv(ra, qs)[:diversity] ≈ r_rma[:diversity]
            @rget r_nma
            @test metadiv(na, qs)[:diversity] ≈ r_nma[:diversity]
            @rget r_rmb
            @test metadiv(rb, qs)[:diversity] ≈ r_rmb[:diversity]
            @rget r_nmb
            @test metadiv(nb, qs)[:diversity] ≈ r_nmb[:diversity]
            @rget r_rmr
            @test metadiv(rr, qs)[:diversity] ≈ r_rmr[:diversity]
            @rget r_nmr
            @test metadiv(nr, qs)[:diversity] ≈ r_nmr[:diversity]
            @rget r_mg
            @test metadiv(g, qs)[:diversity] ≈ r_mg[:diversity]
            
            # Checking out subcommunity diversities
            R"r_rsa = raw_sub_alpha(r_meta, qs)
              r_nsa = norm_sub_alpha(r_meta, qs)
              r_rsb = raw_sub_beta(r_meta, qs)
              r_nsb = norm_sub_beta(r_meta, qs)
              r_rsr = raw_sub_rho(r_meta, qs)
              r_nsr = norm_sub_rho(r_meta, qs)
              r_sg  = sub_gamma(r_meta, qs)"
            @rget r_rsa
            @test subdiv(ra, qs)[:diversity] ≈ r_rsa[:diversity]
            @rget r_nsa
            @test subdiv(na, qs)[:diversity] ≈ r_nsa[:diversity]
            @rget r_rsb
            @test subdiv(rb, qs)[:diversity] ≈ r_rsb[:diversity]
            @rget r_nsb
            @test subdiv(nb, qs)[:diversity] ≈ r_nsb[:diversity]
            @rget r_rsr
            @test subdiv(rr, qs)[:diversity] ≈ r_rsr[:diversity]
            @rget r_nsr
            @test subdiv(nr, qs)[:diversity] ≈ r_nsr[:diversity]
            @rget r_sg
            @test subdiv(g, qs)[:diversity] ≈ r_sg[:diversity]
        end
    end
end

end

end
