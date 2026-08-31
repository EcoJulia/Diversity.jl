# SPDX-License-Identifier: BSD-2-Clause
#
# Cross-validation against other packages: every `test/pkg_*.jl` checks Diversity's results against
# the package it is named for (`test/pkg_Package.jl` validates against `Package`).
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_pkg.jl"])'
#
# **This is the set worth being able to name, and in this package it is not close.**
# `pkg_RCall.jl` reaches `run_rcall.jl`, which is most of the suite's wall-clock on its own: seven
# measures × two scales × random populations, phylogenies and every VCF in `test/data`, each
# iteration crossing the Julia<->R boundary - and, on a cold machine, installing `ape`, `vegan` and
# `rdiversity` from CRAN first. Naming `core_test.jl` instead of this file is the difference between
# a few seconds and several minutes.
#
# **An extra rather than a core set, and the distinction is semantic**: the core sets test *this*
# package against itself, while this one checks it against someone else's - a different question,
# answered against a moving target, and one that only makes sense once the core sets pass. Being an
# `extras_*` also puts it after them in `runtests.jl`, which is where it belongs.
#
# **Not** runnable as a bare script: `Distances`, `RCall` and `StatsBase` are `[extras]` in
# `Project.toml`'s `test` target, and `run_rcall.jl` additionally needs `Phylo` and `PopGen`.
#
# Note: A red cross-validation is a *question*, not automatically a defect here - the reference
# implementation may have changed its convention. See the Gower/vegan note in `run_rcall.jl` and the
# "Before changing what a measure computes" section of `CLAUDE.md` before changing anything.

using Random
using Test
using Diversity
using ParallelTestRunner: find_tests, parse_args, runtests

let pkgbase = map(file -> replace(file, r"pkg_(.*).jl$" => s"\1"),
                  filter(str -> occursin(r"^pkg_.*\.jl$", str),
                         readdir(@__DIR__)))
    if length(pkgbase) > 0
        Random.seed!(1234)
        @info "Cross validation packages:"
        @testset "Cross-validation" begin
            for p in pkgbase
                println("    = $p")
            end
            println()

            runtests(Diversity, parse_args(String[]),
                     testsuite = filter(kv -> startswith(kv.first, "pkg_"),
                                        find_tests(@__DIR__)))
        end
    end
end
