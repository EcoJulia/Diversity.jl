module TestGeneralised
using Compat.Test

using Diversity
using Diversity.ShortNames

metacom = Metacommunity([1 2; 1 2]/6)

# Basic checks for the diversity() function
@testset "diversity()" begin
    @test diversity([metacommunityDiversity], [ρ̄], metacom, 1)[:diversity] ≈ [1.0]
    @test diversity([subcommunityDiversity], [ρ̄], metacom, 1)[:diversity] ≈ [1.0, 1.0]
end

pops = rand(5, 4)
pops = pops / sum(pops)
meta = Metacommunity(pops)
@testset "all-in-one norm/raw_alpha/beta/rho/gamma_div()" begin
    ra = α(meta)
    na = ᾱ(meta)
    rb = β(meta)
    nb = β̄(meta)
    rr = ρ(meta)
    nr = ρ̄(meta)
    g  = Γ(meta)
    # Alpha diversities
    @test raw_sub_alpha(meta, 0)[:diversity] ≈ subdiv(ra, 0)[:diversity]
    @test raw_meta_alpha(meta, 1)[:diversity] ≈ metadiv(ra, 1)[:diversity]
    @test norm_sub_alpha(meta, 2)[:diversity] ≈ subdiv(na, 2)[:diversity]
    @test norm_meta_alpha(meta, 3)[:diversity] ≈ metadiv(na, 3)[:diversity]
    # Beta diversities
    @test raw_sub_beta(meta, 0.0)[:diversity] ≈ subdiv(rb, 0.0)[:diversity]
    @test raw_meta_beta(meta, 1.0)[:diversity] ≈ metadiv(rb, 1.0)[:diversity]
    @test norm_sub_beta(meta, 2.0)[:diversity] ≈ subdiv(nb, 2.0)[:diversity]
    @test norm_meta_beta(meta, 3.0)[:diversity] ≈ metadiv(nb, 3.0)[:diversity]

    @test raw_sub_rho(meta, Inf)[:diversity] ≈ subdiv(rr, Inf)[:diversity]
    @test raw_meta_rho(meta, Inf)[:diversity] ≈ metadiv(rr, Inf)[:diversity]
    @test norm_sub_rho(meta, Inf)[:diversity] ≈ subdiv(nr, Inf)[:diversity]
    @test norm_meta_rho(meta, Inf)[:diversity] ≈ metadiv(nr, Inf)[:diversity]

    # Gamma diversities
    @test sub_gamma(meta, 0.5)[:diversity] ≈ subdiv(g, 0.5)[:diversity]
    @test meta_gamma(meta, 1.5)[:diversity] ≈ metadiv(g, 1.5)[:diversity]

end

end
