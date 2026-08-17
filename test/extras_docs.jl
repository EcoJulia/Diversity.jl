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
# The blocks are extracted and run directly rather than through Documenter, which is why an
# executable block may only `using` packages reachable from the *test* environment. It also means
# this half checks that the code **runs**, not that it produces the output shown: under `@repl`
# Documenter regenerates the output at build time, so there is nothing recorded here to compare
# against. Pinning numbers is `extras_canonical.jl`'s job, not this one's.
#
# **The second half builds the manual with Documenter**, which is a different question and catches
# what the first half structurally cannot: prose. A `[foo](@ref)` pointing at nothing, an
# `@autodocs` that collects no docstring, a page missing from the navigation — none of those are
# code, so none of them can fail above. Before this was added they surfaced only in the
# `Documentation` workflow, after a push.
#
# Both halves are kept, because they fail on different things: Documenter is content with a block
# that emits warnings, and it runs each page's blocks with the manual's own imports rather than
# checking that a page brings its own.

module ExtrasDocs

using Test
using Diversity
using Documenter

# The manual build below runs when this file was asked for by name and when running locally, but not
# on a CI runner that merely reached it as part of the whole suite — the `Documentation` workflow
# already builds the manual, on the same triggers, so doing it again in each of the six test jobs
# checks nothing new. Same "was I asked for?" idiom as `extras_clean.jl`, and for the same reason:
# it needs no environment variable to tell the cases apart. The code blocks above still run
# everywhere, being cheap and platform-dependent in a way the build is not.
const ASKED_FOR = any(a -> occursin("extras_docs", a), ARGS)

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

if !(ASKED_FOR || !haskey(ENV, "RUNNER_OS"))
    @info "Skipping the manual build: this is a CI runner and it was not asked for directly, so " *
          "the `Documentation` workflow builds the manual instead."
else
    @testset "Documentation build" begin
        docsdir = joinpath(@__DIR__, "..", "docs")
        # The same modules, pages and site name `docs/make.jl` publishes with — included rather
        # than repeated, so this cannot drift into checking a different site.
        include(joinpath(docsdir, "config.jl"))
        println()
        @info "Building the manual to check its cross-references ..."

        # `deploydocs` is deliberately not called: this is a check, not a publication. The build
        # goes to a temporary directory so it cannot leave `docs/build` behind for the hygiene
        # tests to find.
        @test isnothing(makedocs(root = docsdir,
                                 modules = DOCS_MODULES,
                                 sitename = DOCS_SITENAME,
                                 pages = DOCS_PAGES,
                                 build = mktempdir()))
    end
end

end
