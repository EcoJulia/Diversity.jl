module ValidateRCall_rdiversity
using Base.Test

using Diversity
using Diversity.ShortNames
using RCall
using DataFrames

# Only run R on unix
if is_unix()
    # Environment variable to avoid boring R package builds
    is_me = haskey(ENV, "SKIP_R_INSTALL") && ENV["SKIP_R_INSTALL"] == "1"
    # Skip the (slow!) R package installation step
    if is_me
        R"library(rdiversity)";
    else
        libdir = mktempdir();
        rcall(Symbol(".libPaths"), libdir);
        rcall(Symbol("install.packages"), "devtools",
              lib=libdir, repos="http://cran.r-project.org");
        @rput libdir
        R"library(devtools, lib.loc=c(libdir, .libPaths()))";
        rcall(Symbol("install_github"), "boydorr/rdiversity", lib=libdir);
        R"library(rdiversity, lib.loc=c(libdir, .libPaths()))";
    end

    # Run diversity comparisons on increasing numbers of types and subcommunities
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
            @testset "Z matrix - $(names[k])" for k in 1:length(Zs)
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
                # Create the metacommunity in R
                r_meta = rcall(:metacommunity, pops, Z)

                # Check out metacommunity diversities
                r_rma = rcall(:raw_meta_alpha, r_meta, qs)
                @test metadiv(ra, qs)[:diversity] ≈ rcopy(r_rma[:diversity])
                r_nma = rcall(:norm_meta_alpha, r_meta, qs)
                @test metadiv(na, qs)[:diversity] ≈ rcopy(r_nma[:diversity])
                r_rmb = rcall(:raw_meta_beta, r_meta, qs)
                @test metadiv(rb, qs)[:diversity] ≈ rcopy(r_rmb[:diversity])
                r_nmb = rcall(:norm_meta_beta, r_meta, qs)
                @test metadiv(nb, qs)[:diversity] ≈ rcopy(r_nmb[:diversity])
                r_rmr = rcall(:raw_meta_rho, r_meta, qs)
                @test metadiv(rr, qs)[:diversity] ≈ rcopy(r_rmr[:diversity])
                r_nmr = rcall(:norm_meta_rho, r_meta, qs)
                @test metadiv(nr, qs)[:diversity] ≈ rcopy(r_nmr[:diversity])
                r_mg = rcall(:meta_gamma, r_meta, qs)
                @test metadiv(g, qs)[:diversity] ≈ rcopy(r_mg[:diversity])
                
                # Check out subcommunity diversities
                r_rsa = rcall(:raw_sub_alpha, r_meta, qs)
                @test subdiv(ra, qs)[:diversity] ≈ rcopy(r_rsa[:diversity])
                r_nsa = rcall(:norm_sub_alpha, r_meta, qs)
                @test subdiv(na, qs)[:diversity] ≈ rcopy(r_nsa[:diversity])
                r_rsb = rcall(:raw_sub_beta, r_meta, qs)
                @test subdiv(rb, qs)[:diversity] ≈ rcopy(r_rsb[:diversity])
                r_nsb = rcall(:norm_sub_beta, r_meta, qs)
                @test subdiv(nb, qs)[:diversity] ≈ rcopy(r_nsb[:diversity])
                r_rsr = rcall(:raw_sub_rho, r_meta, qs)
                @test subdiv(rr, qs)[:diversity] ≈ rcopy(r_rsr[:diversity])
                r_nsr = rcall(:norm_sub_rho, r_meta, qs)
                @test subdiv(nr, qs)[:diversity] ≈ rcopy(r_nsr[:diversity])
                r_sg = rcall(:sub_gamma, r_meta, qs)
                @test subdiv(g, qs)[:diversity] ≈ rcopy(r_sg[:diversity])
            end
        end
    end
end

end
