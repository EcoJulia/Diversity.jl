# SPDX-License-Identifier: BSD-2-Clause

using Test
using Diversity
using ParallelTestRunner: find_tests, parse_args, runtests

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
# The `.jl` is optional, so `test_args = ["extras_clean"]` works too. Any test file may be named -
# `test_Metacommunity.jl` as readily as a whole set.
#
# **The suite is six nameable sets**, which is what lets you run one part rather than all of it:
#
#     core_test  core_ext
#     extras_canonical  extras_clean  extras_docs  extras_pkg
#
# The split is semantic: the **core** sets test this package against itself, the **extras** check
# it against something outside - another package's answers, the blessed results, the documentation,
# the repo's own hygiene.
#
# `extras_pkg` is the one to know about: it cross-validates against R `rdiversity` and `vegan`,
# which is most of the suite's wall-clock and installs CRAN packages on a cold machine. Naming
# `core_test` instead is the difference between seconds and minutes while iterating.
#
# **Running the sets in parallel gives up the ordering guarantee below** - the extras then run even
# when the unit tests are failing, so one broken thing reports as several. If you do, let the first
# invocation get through precompilation before starting the rest, or every process compiles the same
# package at once and they contend.
get!(ENV, "GKSwstype", "100")

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
    corebase = sort(filter(str -> occursin(r"^core_.*\.jl$", str), readdir()))
    extrabase = sort(filter(str -> occursin(r"^extras_.*\.jl$", str),
                            readdir()))

    println()
    @info "Running the core test sets:"
    foreach(f -> println("    = $f"), corebase)
    println()

    @testset "Diversity.jl" begin
        @testset for fn in corebase
            println("    * Running $fn ...")
            include(fn)
        end
    end

    skipextras = get(ENV, "RUNNER_OS", "") == "Windows"
    if skipextras
        println()
        @info "Skipping the extra test suites on a Windows runner."
    elseif !isempty(extrabase)
        # Three groups, in this order: the serial extras, then the concurrent ones, then hygiene.
        #
        # `extras_clean` being **last of all** rather than beside its siblings is load-bearing, not
        # tidiness. It fails on any unstaged change to a tracked file - the normal state of a
        # working tree mid-task - and on a stale `dateModified` in `codemeta.json`, which goes stale
        # overnight. A `@testset` throws at its end, so while it sat in the serial group that
        # routine failure aborted the run before `extras_docs` had even started, and the
        # documentation went unchecked locally. Moving it within the serial group would not have
        # helped: the throw comes from the enclosing testset, whichever member failed.
        #
        # The cost of the swap is the mirror case - a failing `extras_docs` now stops `extras_clean`
        # from running - which is much the better trade, since a broken docs build is a real defect
        # while a dirty tree is not.
        parallelextras = ["extras_docs", "extras_examples", "extras_notebooks"]
        lastextras = ["extras_clean"]
        setname(fn) = chop(fn, tail = 3)
        serialextras = filter(fn -> setname(fn) ∉ parallelextras &&
                                    setname(fn) ∉ lastextras, extrabase)
        finalextras = filter(fn -> setname(fn) ∈ lastextras, extrabase)

        println()
        @info "Running the extra test suites, in this order:"
        foreach(f -> println("    = $f"), serialextras)
        foreach(f -> println("    = $f (concurrently)"),
                filter(fn -> setname(fn) ∈ parallelextras, extrabase))
        foreach(f -> println("    = $f (last)"), finalextras)
        println()

        if !isempty(serialextras)
            @testset "Serial extras..." begin
                @testset for fn in serialextras
                    println("    * Running $fn ...")
                    include(fn)
                end
            end
        end

        suite = filter(kv -> kv.first in parallelextras, find_tests(@__DIR__))
        if !isempty(suite)
            println()
            @info "Running these concurrently: " *
                  join(sort(collect(keys(suite))), ", ")
            println()
            runtests(Diversity, parse_args(String[]), testsuite = suite)
        end

        if !isempty(finalextras)
            println()
            @testset "Hygiene, last of all..." begin
                @testset for fn in finalextras
                    println("    * Running $fn ...")
                    include(fn)
                end
            end
        end
    end
end
