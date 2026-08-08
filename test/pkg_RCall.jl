# SPDX-License-Identifier: BSD-2-Clause

module ValidateRCall

using Test

# Environment variable to avoid boring R package builds
mustCrossvalidate = haskey(ENV, "JULIA_MUST_CROSSVALIDATE") &&
                    ENV["JULIA_MUST_CROSSVALIDATE"] == "1"

# Only run R on unix or when R is installed because JULIA_MUST_CROSSVALIDATE is set to 1
global skipR = !mustCrossvalidate && !Sys.isunix()
try
    skipR && error("Skipping R testing...")
    using RCall
    global skipR = false
catch
    global skipR = true
    @warn "R or appropriate Phylo package not installed, skipping R cross-validation."
end

if skipR
    # ⚠️ Say so in the test summary, not only in the log. Without this the file contributes *no* test
    # results at all, so a run that skipped the entire R cross-validation — nine thousand assertions
    # — looks exactly like one that passed it.
    @info "Skipping R cross-validation. Set JULIA_MUST_CROSSVALIDATE=1 to make this an error " *
          "rather than a skip."
    @testset "R cross-validation" begin
        @test_broken "R cross-validation was skipped" == ""
    end
else
    include("run_rcall.jl")
end

end
