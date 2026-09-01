# SPDX-License-Identifier: BSD-2-Clause

module CleanJuliaFormatter
using Test
using Diversity
using Git
using JuliaFormatter

include("GitUtils.jl")
using .GitUtils

@testset "JuliaFormatter" begin
    git_dir = readchomp(`$(Git.git()) rev-parse --show-toplevel`)
    @test_nowarn format(Diversity)
    @test is_repo_clean(git_dir, strict = haskey(ENV, "RUNNER_OS"))
end

end
