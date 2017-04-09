module ValidateRCall_rdiversity
using Base.Test

using Diversity
using Diversity.ShortNames
using RCall
using DataFrames

R"install.packages('devtools', repos='http://cran.r-project.org')
  devtools::install_github('boydorr/rdiversity')
  library(rdiversity)"
        
@testset "RCall - rdiversity" begin
    @testset "Random rdiversity $i" for i in 1:20
        types = rand(1:(i*10))
        sc = rand(1:(i*10))
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
        pops /= sum(pops)
        meta_sym = Metacommunity(pops, Z)
        meta_asym = Metacommunity(pops, Zasym)
        qs = sort([rand(7)*10..., 0, 1, Inf])
        sra = α(meta_sym)
        sna = ᾱ(meta_sym)
        srb = β(meta_sym)
        snb = β̄(meta_sym)
        srr = ρ(meta_sym)
        snr = ρ̄(meta_sym)
        sg  = Γ(meta_sym)
        ara = α(meta_asym)
        ana = ᾱ(meta_asym)
        arb = β(meta_asym)
        anb = β̄(meta_asym)
        arr = ρ(meta_asym)
        anr = ρ̄(meta_asym)
        ag  = Γ(meta_asym)
        @rput pops
        @rput qs
        @rput Z
        @rput Zasym
        # Checking out metacommunity diversities
        R"meta_sym = metacommunity(pops, Z)
          meta_asym = metacommunity(pops, Zasym)
          r_srma = raw_meta_alpha(meta_sym, qs)
          r_snma = norm_meta_alpha(meta_sym, qs)
          r_srmb = raw_meta_beta(meta_sym, qs)
          r_snmb = norm_meta_beta(meta_sym, qs)
          r_srmr = raw_meta_rho(meta_sym, qs)
          r_snmr = norm_meta_rho(meta_sym, qs)
          r_smg  = meta_gamma(meta_sym, qs)
          r_arma = raw_meta_alpha(meta_asym, qs)
          r_anma = norm_meta_alpha(meta_asym, qs)
          r_armb = raw_meta_beta(meta_asym, qs)
          r_anmb = norm_meta_beta(meta_asym, qs)
          r_armr = raw_meta_rho(meta_asym, qs)
          r_anmr = norm_meta_rho(meta_asym, qs)
          r_amg  = meta_gamma(meta_asym, qs)"
        @rget r_srma
        @test metadiv(sra, qs)[:diversity] ≈ r_srma[:diversity]
        @rget r_snma
        @test metadiv(sna, qs)[:diversity] ≈ r_snma[:diversity]
        @rget r_arma
        @test metadiv(ara, qs)[:diversity] ≈ r_arma[:diversity]
        @rget r_anma
        @test metadiv(ana, qs)[:diversity] ≈ r_anma[:diversity]
        @rget r_srmb
        @test metadiv(srb, qs)[:diversity] ≈ r_srmb[:diversity]
        @rget r_snmb
        @test metadiv(snb, qs)[:diversity] ≈ r_snmb[:diversity]
        @rget r_armb
        @test metadiv(arb, qs)[:diversity] ≈ r_armb[:diversity]
        @rget r_anmb
        @test metadiv(anb, qs)[:diversity] ≈ r_anmb[:diversity]
        @rget r_srmr
        @test metadiv(srr, qs)[:diversity] ≈ r_srmr[:diversity]
        @rget r_snmr
        @test metadiv(snr, qs)[:diversity] ≈ r_snmr[:diversity]
        @rget r_armr
        @test metadiv(arr, qs)[:diversity] ≈ r_armr[:diversity]
        @rget r_anmr
        @test metadiv(anr, qs)[:diversity] ≈ r_anmr[:diversity]
        @rget r_smg
        @test metadiv(sg, qs)[:diversity] ≈ r_smg[:diversity]
        @rget r_amg
        @test metadiv(ag, qs)[:diversity] ≈ r_amg[:diversity]

        # Checking out subcommunity diversities
        R"r_srsa = raw_sub_alpha(meta_sym, qs)
          r_snsa = norm_sub_alpha(meta_sym, qs)
          r_srsb = raw_sub_beta(meta_sym, qs)
          r_snsb = norm_sub_beta(meta_sym, qs)
          r_srsr = raw_sub_rho(meta_sym, qs)
          r_snsr = norm_sub_rho(meta_sym, qs)
          r_ssg  = sub_gamma(meta_sym, qs)
          r_arsa = raw_sub_alpha(meta_asym, qs)
          r_ansa = norm_sub_alpha(meta_asym, qs)
          r_arsb = raw_sub_beta(meta_asym, qs)
          r_ansb = norm_sub_beta(meta_asym, qs)
          r_arsr = raw_sub_rho(meta_asym, qs)
          r_ansr = norm_sub_rho(meta_asym, qs)
          r_asg  = sub_gamma(meta_asym, qs)"
        @rget r_srsa
        @test subdiv(sra, qs)[:diversity] ≈ r_srsa[:diversity]
        @rget r_snsa
        @test subdiv(sna, qs)[:diversity] ≈ r_snsa[:diversity]
        @rget r_arsa
        @test subdiv(ara, qs)[:diversity] ≈ r_arsa[:diversity]
        @rget r_ansa
        @test subdiv(ana, qs)[:diversity] ≈ r_ansa[:diversity]
        @rget r_srsb
        @test subdiv(srb, qs)[:diversity] ≈ r_srsb[:diversity]
        @rget r_snsb
        @test subdiv(snb, qs)[:diversity] ≈ r_snsb[:diversity]
        @rget r_arsb
        @test subdiv(arb, qs)[:diversity] ≈ r_arsb[:diversity]
        @rget r_ansb
        @test subdiv(anb, qs)[:diversity] ≈ r_ansb[:diversity]
        @rget r_srsr
        @test subdiv(srr, qs)[:diversity] ≈ r_srsr[:diversity]
        @rget r_snsr
        @test subdiv(snr, qs)[:diversity] ≈ r_snsr[:diversity]
        @rget r_arsr
        @test subdiv(arr, qs)[:diversity] ≈ r_arsr[:diversity]
        @rget r_ansr
        @test subdiv(anr, qs)[:diversity] ≈ r_ansr[:diversity]
        @rget r_ssg
        @test subdiv(sg, qs)[:diversity] ≈ r_ssg[:diversity]
        @rget r_asg
        @test subdiv(ag, qs)[:diversity] ≈ r_asg[:diversity]

    end
end

end
