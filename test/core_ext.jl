# SPDX-License-Identifier: BSD-2-Clause
#
# The extension tests: every `test/ext_*.jl`, which test the matching `ext/*.jl`.
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["core_ext.jl"])'
#
# **Not** runnable as a bare script, and that is a dependency fact rather than an oversight: every
# trigger package here - `AxisArrays`, `BioSequences`, `Phylo`, `PopGen` - is a `[weakdeps]` of the
# package and an `[extras]` of `Project.toml`'s `test` target, so only `Pkg.test` provisions them.
# Loading each is what activates the extension under test in the first place.

using Random
using Test
using Diversity

# Identify files in test/ that are testing matching files in ext/
#  - ext/SourceExt.jl will be matched by test/ext_SourceExt.jl
let filebase = String[]
    for (root, dirs, files) in walkdir(joinpath(@__DIR__, "..", "ext"))
        append!(filebase,
                map(file -> replace(file, r"(.*).jl" => s"\1"),
                    filter(file -> occursin(r".*\.jl", file), files)))
    end

    extbase = map(file -> replace(file, r"ext_(.*).jl" => s"\1"),
                  filter(str -> occursin(r"^ext_.*\.jl$", str),
                         readdir(@__DIR__)))

    # Identify tests with no matching file
    superfluous = filter(f -> f ∉ filebase, extbase)
    if length(superfluous) > 0
        println()
        @info "Potentially superfluous extension tests:"
        for f in superfluous
            println("    + $f.jl")
        end
        println()
    end

    # Identify files with no matching test
    notest = filter(f -> f ∉ extbase, filebase)
    if length(notest) > 0
        println()
        @info "Potentially missing extension tests:"
        for f in notest
            println("    - $f.jl")
        end
        println()
    end

    Random.seed!(1234)

    @testset "Extension tests" begin
        println()
        @info "Running tests for extensions:"
        for t in extbase
            println("    = $t.jl")
        end
        println()

        @info "Running extension tests..."
        @testset for t in extbase
            fn = "ext_$t.jl"
            println("    * Testing $t.jl extension...")
            include(joinpath(@__DIR__, fn))
        end
    end
end
