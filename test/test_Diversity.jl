module TestDiversity
using Base.Test

using Diversity

@testset "Deprecations" begin
    # Only run deprecation checks on macs
    if is_apple()
        @test_warn "deprecated" subcommunityalpha([1:3...]/6, 1)
        @test_warn "deprecated" subcommunityalphabar([1:3...]/6, 1)
        @test_warn "deprecated" subcommunitybeta([1:3...]/6, 1)
        @test_warn "deprecated" subcommunitybetabar([1:3...]/6, 1)
        @test_warn "deprecated" subcommunityrho([1:3...]/6, 1)
        @test_warn "deprecated" subcommunityrhobar([1:3...]/6, 1)
        @test_warn "deprecated" subcommunitygamma([1:3...]/6, 1)
        @test_warn "deprecated" subcommunitygammabar([1:3...]/6, 1)

        @test_warn "deprecated" supercommunityA([1:3...]/6, 1)
        @test_warn "deprecated" supercommunityAbar([1:3...]/6, 1)
        @test_warn "deprecated" supercommunityB([1:3...]/6, 1)
        @test_warn "deprecated" supercommunityBbar([1:3...]/6, 1)
        @test_warn "deprecated" supercommunityR([1:3...]/6, 1)
        @test_warn "deprecated" supercommunityRbar([1:3...]/6, 1)
        @test_warn "deprecated" supercommunityG([1:3...]/6, 1)
        @test_warn "deprecated" supercommunityGbar([1:3...]/6, 1)

        @test_warn "deprecated" ecosystemA([1:3...]/6, 1)
        @test_warn "deprecated" ecosystemAbar([1:3...]/6, 1)
        @test_warn "deprecated" ecosystemB([1:3...]/6, 1)
        @test_warn "deprecated" ecosystemBbar([1:3...]/6, 1)
        @test_warn "deprecated" ecosystemR([1:3...]/6, 1)
        @test_warn "deprecated" ecosystemRbar([1:3...]/6, 1)
        @test_warn "deprecated" ecosystemG([1:3...]/6, 1)
        @test_warn "deprecated" ecosystemGbar([1:3...]/6, 1)
    end
end

end
