# SPDX-License-Identifier: BSD-2-Clause
#
# Run the executable code in `docs/src/*.md`, so a documentation page cannot go on showing code that
# no longer works.
#
# On its own:
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_docs.jl"])'
#
# **Not** runnable as a bare script, and that is a dependency fact rather than an oversight:
# `phylogenetics.md` and `genetics.md` `using` `Phylo`, `PopGen` and `BioSequences`, which are
# `[extras]` in `Project.toml`'s `test` target.
#
# **A ```@repl fence is the single source of truth for "this code is checked".** Documenter runs
# those blocks when it builds the site; this file runs the same blocks in the test suite. A plain
# ```julia fence stays illustrative and is deliberately *not* run. One marker, two runners, and no
# second list to keep in step.
#
# This does **not** call Documenter, and cannot: `Documenter` lives in `docs/Project.toml`, not in
# `Project.toml`'s `test` target, so `makedocs` is unavailable here. The blocks are extracted and run
# directly instead — which also means an executable block may only `using` packages reachable from
# the *test* environment.
#
# It also means this file checks that the code **runs**, not that it produces the output shown:
# under `@repl` Documenter regenerates the output at build time, so there is nothing recorded here to
# compare against. Pinning numbers is `extras_canonical.jl`'s job, not this one's.

module ExtrasDocs

using Test
using Diversity

# Pages that plot need GR told there is no display, exactly as `docs/make.jl` does — otherwise a
# headless runner fails on the first figure.
get!(ENV, "GKSwstype", "100")

# The fence languages Documenter *executes* while building a page. `@example` and `@setup` run as
# scripts (the latter without showing its code); `@repl` runs line by line and shows a prompt. The
# rest of Documenter's blocks — `@meta`, `@docs`, `@autodocs`, `@index`, `@contents`, `@raw` — are
# directives rather than code, and `jldoctest` is checked against its own recorded output by
# Documenter itself, so running it here would duplicate that badly (we would execute it but compare
# nothing).
const EXECUTED = ("@example", "@setup", "@repl")

# One page's executable code, split into sandboxes **exactly as Documenter groups them**: blocks
# sharing a name share a module and run in order, while an anonymous block gets a module to itself.
# Matching that rule is what makes this runner and the docs build agree by construction — a page that
# passes here for the wrong reason (state leaking between blocks that Documenter keeps apart) would
# fail there.
#
# Returns name => concatenated code, in first-appearance order.
function _sandboxes(path::AbstractString)
    order = String[]
    code = Dict{String, Vector{String}}()
    anonymous = 0
    infence = false
    key = ""
    buffer = String[]
    for line in eachline(path)
        if infence
            if startswith(line, "```")
                infence = false
                push!(code[key], join(buffer, "\n"))
                empty!(buffer)
            else
                push!(buffer, line)
            end
        elseif startswith(line, "```")
            # ```@repl name — the language is the first word, the sandbox name the rest.
            spec = split(strip(chopprefix(line, "```")); limit = 2)
            isempty(spec) && continue
            first(spec) in EXECUTED || continue
            name = length(spec) == 2 ? strip(spec[2]) : ""
            if isempty(name)
                anonymous += 1
                name = "anonymous-$anonymous"
            end
            key = name
            haskey(code, key) || (push!(order, key); code[key] = String[])
            infence = true
        end
    end
    # An unterminated fence means the page is malformed; say so rather than silently running a
    # truncated block or none at all.
    infence && error("unterminated code fence in $(basename(path))")
    return [name => join(code[name], "\n") for name in order]
end

# Run one sandbox and assert two separate things: that it does not throw, and that it produces no
# warning or error on stderr.
#
# It deliberately does **not** use `@test_nowarn`, which fails on *any* stderr output including
# `@info`. That is too strict for documentation: a page is entitled to call a package that logs
# — `SpatialEcology` announces "Matrix data assumed to be presence-absence" whenever an assemblage is
# built — and forbidding that would mean either hiding the call or dropping the example. A *warning*
# still fails, because a documentation example that warns is usually a documentation example doing
# something wrong.
function _runblock(sandbox, source, label)
    ok, log = mktemp() do path, io
        result = redirect_stderr(io) do
            try
                include_string(sandbox, source, label)
                true
            catch e
                @error "documentation block failed" label exception = e
                false
            end
        end
        flush(io)
        return result, read(path, String)
    end
    @test ok
    @test !occursin("Warning:", log) && !occursin("Error:", log)
    return nothing
end

@testset "Documentation code" begin
    docsdir = joinpath(@__DIR__, "..", "docs", "src")
    pages = sort(filter(f -> endswith(f, ".md"), readdir(docsdir)))
    println()
    @info "Running the executable code in docs/src ..."
    total = 0
    for page in pages
        sandboxes = _sandboxes(joinpath(docsdir, page))
        isempty(sandboxes) && continue
        total += length(sandboxes)
        println("    * $page — $(length(sandboxes)) executable block group(s) ...")
        @testset "$page" begin
            for (name, source) in sandboxes
                # A fresh, bare module per sandbox: the page must bring its own `using` statements,
                # which is the point — a page whose imports only work because the test suite had
                # already loaded something is a page a reader cannot follow.
                sandbox = Module(Symbol("Docs_", replace(page, r"\W" => "_"),
                                        "_",
                                        replace(name, r"\W" => "_")))
                @testset "$name" begin
                    _runblock(sandbox, source, "$page [$name]")
                end
            end
        end
    end
    # The check that this file is doing anything at all. A regex that quietly matches nothing
    # reports success just as loudly as one that works, and this suite exists precisely because
    # unexecuted documentation rots invisibly — so a run that executed no code is a failure.
    @test total > 0
end

end
