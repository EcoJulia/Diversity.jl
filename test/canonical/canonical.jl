# SPDX-License-Identifier: BSD-2-Clause
#
# The blessing machinery shared by every canonical test. See README.md in this directory for what a
# canonical test is for and how to re-bless one.

module Canonical

using Test
using TOML

# Note: The function is `blessed`, not `canonical`, and that is not a style choice: **`BioSequences`
# exports `canonical`** (the canonical orientation of a k-mer), so a canonical test file that loads it
# — `test_types.jl` does — gets an ambiguity rather than either function, reported as a bare
# `UndefVarError`. Checked against every package the test target can load; `blessed` is free in all of
# them. Do not "restore parity" with EcoSISTEM by renaming this back.
export blessed, canonical_reference, writereference, blessing

# Where the blessed numbers live. One file for the whole directory, so a re-blessing produces a single
# reviewable diff rather than a scatter of them.
const REFERENCE = joinpath(@__DIR__, "reference.toml")

# Note: Blessed values are written through **on every call**, not accumulated and flushed at the end.
# That looks wasteful and is deliberate: each canonical test file `include`s this file into its own
# module, so every file gets its *own* copy of any in-memory state — an accumulate-then-flush design
# silently blesses nothing at all, because the runner writes out its own empty dict rather than the
# values the test files recorded into theirs. Writing through cannot have that bug.
const RECORDED = Dict{String, Any}()

"""
    blessing()

Is this a re-blessing run? True when `DIVERSITY_BLESS=true`, in which case the canonical tests record
their results instead of checking them.
"""
blessing() = get(ENV, "DIVERSITY_BLESS", "false") == "true"

# The blessed values, or an empty dict the first time. Read once and cached, since every call needs it.
const _CACHE = Ref{Union{Nothing, Dict{String, Any}}}(nothing)
function canonical_reference()
    isnothing(_CACHE[]) &&
        (_CACHE[] = isfile(REFERENCE) ? TOML.parsefile(REFERENCE) :
                    Dict{String, Any}())
    return _CACHE[]
end

# TOML holds numbers and flat arrays of them. `TOML.print` **errors outright on a `Matrix`**, so a
# similarity matrix has to be flattened — and refusing it here, with the fix in the message, is much
# clearer than letting the write fail deep inside a re-blessing run. Flattening at the *call site* is
# also the point: it leaves the shape asserted in the test, where a reader can see it, rather than
# implied by the length of an array in the reference file.
function _plain(name, value)
    value isa Real && return float(value)
    value isa AbstractVector{<:Real} && return float.(collect(value))
    value isa AbstractArray{<:Real} &&
        return error("canonical value `$name` is a $(ndims(value))-dimensional array; TOML stores " *
                     "flat arrays only. Flatten it explicitly — `vec(Z)` — and assert its shape " *
                     "separately in the test.")
    return error("canonical value `$name` is a $(typeof(value)); a blessed value must be a real " *
                 "number or a vector of them.")
end

"""
    blessed(name, value; rtol = 1e-8)

Compare `value` against the blessed result for `name`, or record it when re-blessing.

Flatten matrices before calling — `vec(Z)` — and assert the shape in the test instead.

`rtol` is deliberately tight by default. A canonical test exists to notice change, so a loose
tolerance defeats it; widen it only where a result is genuinely only reproducible to fewer digits, and
say why at the call site.
"""
function blessed(name::AbstractString, value; rtol = 1e-8)
    key = String(name)
    plain = _plain(key, value)
    RECORDED[key] = plain
    if blessing()
        _writethrough(key, plain)
        return @test true                    # nothing to compare against; this run defines it
    end
    ref = canonical_reference()
    if !haskey(ref, key)
        return @test_broken "no blessed value for `$key` — run the canonical suite with " *
                            "DIVERSITY_BLESS=true to record one" == ""
    end
    return @test isapprox(plain, ref[key]; rtol = rtol)
end

# Merge one blessed value into the reference file. Read-modify-write per call, for the reason given
# at `RECORDED` above; the file is small and blessing is rare.
function _writethrough(key, plain)
    merged = merge(canonical_reference(), Dict(key => plain))
    _CACHE[] = merged
    _write(merged)
    return nothing
end

function _write(merged)
    open(REFERENCE, "w") do io
        println(io,
                "# Blessed canonical results — regenerate with DIVERSITY_BLESS=true, and read")
        println(io,
                "# test/canonical/README.md before committing a change to this file.")
        return TOML.print(io, merged; sorted = true)
    end
    return nothing
end

"""
    writereference()

Write everything recorded this run to `reference.toml`. Call once, after all canonical tests.

Note: **Merges rather than replaces.** A run that executed only some of the canonical files would
otherwise silently delete the blessed values of the rest, turning a partial re-blessing into a
wholesale loss — the sort of damage that only shows up much later, as a test that stopped checking
anything.
"""
function writereference()
    blessing() || return nothing
    merged = merge(canonical_reference(), RECORDED)
    _write(merged)
    @info "blessed $(length(RECORDED)) canonical value(s) into $(REFERENCE)"
    return nothing
end

end
