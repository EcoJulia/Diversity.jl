# SPDX-License-Identifier: BSD-2-Clause

module CleanRSMD
using Test
using Git
using Logging
using ResearchSoftwareMetadata

include("GitUtils.jl")
using .GitUtils

@testset "RSMD" begin
    git_dir = readchomp(`$(Git.git()) rev-parse --show-toplevel`)
    @test isnothing(ResearchSoftwareMetadata.crosswalk())
    global_logger(SimpleLogger(stderr, Logging.Warn))
    @test_nowarn ResearchSoftwareMetadata.crosswalk()
    global_logger(SimpleLogger(stderr, Logging.Info))
    @test is_repo_clean(git_dir, strict = haskey(ENV, "RUNNER_OS"))
end

end
