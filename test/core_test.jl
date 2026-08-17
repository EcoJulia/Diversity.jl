# SPDX-License-Identifier: BSD-2-Clause
#
# The unit tests: every `test/test_*.jl`, which test the matching `src/*.jl`.
#
# Run this set on its own — everything the package computes, with none of the cross-validation:
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["core_test.jl"])'
#
# That is the point of this file existing. Any *single* file could already be named
# (`test_args = ["test_Metacommunity.jl"]`), but there was no way to say "all the unit tests and
# nothing else" short of naming all of them — so the only alternative to one file was the whole
# suite, `extras_pkg.jl` and its R cross-validation included.
#
# Most of these files also run as bare scripts (`cd test && julia --project=.. test_Types.jl`),
# because they need only the package's own dependencies. The exception is `test_EcoBase.jl`, which
# needs `SpatialEcology` and `CSV` — `[extras]` in `Project.toml`'s `test` target, and so reachable
# only through `Pkg.test`.

using Random
using Test
using Diversity
using ParallelTestRunner: find_tests, parse_args, runtests

# Identify files in test/ that are testing matching files in src/
#  - src/Source.jl will be matched by test/test_Source.jl
let filebase = String[]
    for (root, dirs, files) in walkdir(joinpath(@__DIR__, "..", "src"))
        append!(filebase,
                map(file -> replace(file, r"(.*).jl" => s"\1"),
                    filter(file -> occursin(r".*\.jl", file), files)))
    end

    testbase = map(file -> replace(file, r"test_(.*).jl" => s"\1"),
                   filter(str -> occursin(r"^test_.*\.jl$", str),
                          readdir(@__DIR__)))

    # Identify tests with no matching file
    superfluous = filter(f -> f ∉ filebase, testbase)
    if length(superfluous) > 0
        println()
        @info "Potentially superfluous tests:"
        for f in superfluous
            println("    + $f.jl")
        end
        println()
    end

    # Identify files with no matching test
    notest = filter(f -> f ∉ testbase, filebase)
    if length(notest) > 0
        println()
        @info "Potentially missing tests:"
        for f in notest
            println("    - $f.jl")
        end
        println()
    end

    # Seeded here rather than only in `runtests.jl`, or this file run on its own would not be
    # reproducible — which is most of why it exists.
    Random.seed!(1234)

    @testset "Unit tests" begin
        @test isfile(pkgdir(Diversity, "test", "runtests.jl"))
        println()
        @info "Running tests for files:"
        for t in testbase
            println("    = $t.jl")
        end
        println()

        runtests(Diversity, parse_args(String[]),
                 testsuite = filter(kv -> startswith(kv.first, "test_"),
                                    find_tests(@__DIR__)))
    end
end
