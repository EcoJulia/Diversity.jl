# SPDX-License-Identifier: BSD-2-Clause
#
# Run the canonical suite — the blessed-result tests in `test/canonical/`.
#
# On its own:
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_canonical.jl"])'
#
# To re-bless after an intended change (read `test/canonical/README.md` first):
#
#     DIVERSITY_BLESS=true julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_canonical.jl"])'
#
# **Why this set exists at all, given the R cross-validation.** `pkg_RCall.jl` is a far stronger
# check — it compares against an independent implementation of the same specification — but it only
# runs where R is installed, which `testing.yaml` arranges on macOS alone. On Ubuntu and Windows
# nothing otherwise notices that a number this package computes has changed. These blessed values are
# cheap, need no R, and run everywhere.
#
# **Not** runnable as a bare script: the type-construction file needs `Phylo`, `PopGen` and
# `BioSequences`, all `[extras]` in `Project.toml`'s `test` target.
#
# Discovery here is by glob (`canonical/test_*.jl`), and that is a trap worth knowing: a file
# dropped into `canonical/` under any other name is never run, and will rot silently.

module ExtrasCanonical

using Test

include(joinpath(@__DIR__, "canonical", "canonical.jl"))
using .Canonical

@testset "Canonical results" begin
    dir = joinpath(@__DIR__, "canonical")
    files = sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"),
                        readdir(dir)))
    println()
    @info "Running canonical tests" * (blessing() ? " — RE-BLESSING" : "")
    for f in files
        println("    * ", f, " ...")
        include(joinpath(dir, f))
    end
    # Written once, after every file, and merged rather than replaced — see `writereference`.
    writereference()
end

end
