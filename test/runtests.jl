# SPDX-License-Identifier: BSD-2-Clause

using Test
using Diversity

# A test argument names one test file to run *instead of* the whole suite:
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_clean.jl"])'
#
# Going through `Pkg.test` rather than running the file directly is the whole point: it is what
# provisions the test environment. `Git`, `JuliaFormatter`, `ResearchSoftwareMetadata`, `RCall`,
# `Phylo`, `PopGen` and the rest are `[extras]` in `Project.toml`'s `test` target and nothing else
# supplies them, so a bare `julia test/extras_clean.jl` dies on `using Git`. `Pkg.test` gets it
# right by construction; anything else reconstructs it and drifts.
#
# The `.jl` is optional, so `test_args = ["extras_clean"]` works too. Any test file may be named —
# `test_Metacommunity.jl` as readily as a whole set.
#
# **The suite is six nameable sets**, which is what lets you run one part rather than all of it:
#
#     core_test  core_ext
#     extras_canonical  extras_clean  extras_docs  extras_pkg
#
# The split is semantic: the **core** sets test this package against itself, the **extras** check
# it against something outside — another package's answers, the blessed results, the documentation,
# the repo's own hygiene.
#
# `extras_pkg` is the one to know about: it cross-validates against R `rdiversity` and `vegan`,
# which is most of the suite's wall-clock and installs CRAN packages on a cold machine. Naming
# `core_test` instead is the difference between seconds and minutes while iterating.
#
# **Running the sets in parallel gives up the ordering guarantee below** — the extras then run even
# when the unit tests are failing, so one broken thing reports as several. If you do, let the first
# invocation get through precompilation before starting the rest, or every process compiles the same
# package at once and they contend.

requested = map(a -> endswith(a, ".jl") ? a : a * ".jl", ARGS)
for fn in requested
    isfile(joinpath(@__DIR__, fn)) ||
        error("`$fn` was asked for with `test_args`, but there is no such file in `test/`.")
end

if !isempty(requested)
    @info "Running only the requested test file(s): " * join(requested, ", ")
    @testset for fn in requested
        println("    * Running $fn ...")
        include(fn)
    end
else
    # Two loops, and nothing else. Each `core_*.jl` and `extras_*.jl` is a standalone set that can be
    # run on its own by name (see above); this file only decides the order they go in.
    #
    # The extras run **after** the core sets deliberately: a failing `@testset` throws at its end,
    # so the extras are reached only once the unit and extension tests pass. There is no point
    # cross-validating a broken package against R, or blessing results it computed wrongly.
    #
    corebase = sort(filter(str -> occursin(r"^core_.*\.jl$", str),
                           readdir(@__DIR__)))
    extrabase = sort(filter(str -> occursin(r"^extras_.*\.jl$", str),
                            readdir(@__DIR__)))

    println()
    @info "Running the core test sets:"
    foreach(f -> println("    = $f"), corebase)
    println()

    @testset "Diversity.jl" begin
        @testset for fn in corebase
            println("    * Running $fn ...")
            include(joinpath(@__DIR__, fn))
        end
    end

    if !isempty(extrabase)
        println()
        @info "Running the extra test suites:"
        foreach(f -> println("    = $f"), extrabase)
        println()

        # Wrapped in an enclosing testset, exactly as the core loop is, and it is **not**
        # decoration: a failing `@testset` throws when it is the *outermost* one, so a bare
        # `@testset for` here would abort the loop at the first set that failed.
        @testset "Extras" begin
            @testset for fn in extrabase
                println("    * Running $fn ...")
                include(joinpath(@__DIR__, fn))
            end
        end
    end
end
