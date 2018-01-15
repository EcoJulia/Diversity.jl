module TestDiversity
using Compat.Test

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

        @test_warn "deprecated" rawmetaalpha([1:3...]/6, 1)
        @test_warn "deprecated" normmetaalpha([1:3...]/6, 1)
        @test_warn "deprecated" rawmetabeta([1:3...]/6, 1)
        @test_warn "deprecated" normmetabeta([1:3...]/6, 1)
        @test_warn "deprecated" rawmetarho([1:3...]/6, 1)
        @test_warn "deprecated" normmetarho([1:3...]/6, 1)
        @test_warn "deprecated" metagamma([1:3...]/6, 1)
    end
end

end
