# Canonical tests

A canonical test records what the package **actually produced** and checks that it keeps producing it.
It answers a different question from the rest of the suite: a unit test asks *is this right?*, a
canonical test asks *has this changed?* Both matter, and neither substitutes for the other.

Run them:

```julia
julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_canonical.jl"])'
```

## Why these exist, given `pkg_RCall.jl`

The R cross-validation is a strictly stronger check — it compares against
[rdiversity](https://github.com/boydorr/rdiversity), an independent implementation of the same
specification, rather than against this package's own past output.

But it only runs where R is installed, and `testing.yaml` arranges that on **macOS alone**. On Ubuntu
and Windows, nothing else notices that a number this package computes has changed. These blessed
values need no R, run in a second, and run everywhere. They are the floor, not the ceiling.

## Re-blessing

When you change what the package computes and the change is *intended*, the blessed numbers must be
updated:

```julia
DIVERSITY_BLESS=true julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_canonical.jl"])'
```

This rewrites `reference.toml`. **The diff to that file is the deliverable** — it is the
machine-checked statement of what your change did to the output, and it belongs in the pull request
alongside the code. Review it before committing:

- Did exactly the results you expected to move, move?
- Did anything move that you did not expect? That is the finding, not the noise.
- Are the new numbers still sensible — positive, finite, non-increasing in `q`?

⚠️ **Re-blessing is not a way to make a failing test pass.** A canonical failure means the output
changed; your job is to explain *why* before recording the new value. If you cannot say why, do not
bless it.

🔴 And if a change moves these numbers, it will very likely move rdiversity's too — so it is a change
to the *specification*, not just to this package. Raise it rather than blessing past it.

## The two kinds, and why there are two

| file | input | what only it can catch |
|---|---|---|
| `test_measures.jl` | a fixed abundance matrix and a fixed `Z` | a change in the **measures**: the ordinariness ratios, the power-mean orders, the raw/normalised weighting |
| `test_types.jl` | a fixed tree, the shipped VCF, fixed sequences | a change in **type construction**: how a `PhyloBranches`, `GeneticVCF` or `GeneticFASTA` turns its input into a similarity matrix before any measure sees it |

Keeping them apart is the point: when a number moves, which file moved tells you immediately whether
to look at `src/DiversityMeasure.jl` or at an extension.

## 🔴 `test_paper.jl` is a third kind, and the strongest

Everything above records what this package produced. `test_paper.jl` records what the **paper says it
must produce** — the communities its appendices work right through, with every measure stated. So it
uses **plain `@test`s and no blessed values at all**, and `DIVERSITY_BLESS=true` cannot touch it.

That is the point. A blessed value asks *has this changed?* and can be silenced by re-blessing;
`test_paper.jl` asks *is this still the published framework?* and cannot. If it fails, either the
package has stopped implementing the specification or the specification has moved — and both are
findings, not something to record and move past.

⚠️ Add to it rather than blessing around it whenever a result has a published value to check against.
And name the examples rather than numbering them: supplementary section numbers move between versions
of a paper.

⭐ `test_types.jl` is also the **only numeric gate on the extensions**. `ext_*.jl` checks them against
hand-computed expectations for one tiny input; nothing else records what they produce.

## Writing one

```julia
include("canonical.jl")
using .Canonical

blessed("measures/unique/Gamma/meta_q0", metadiv(Γ(meta), 0)[1, :diversity])
```

- **Flatten matrices at the call site** — `vec(Z)`. `TOML.print` errors outright on a `Matrix`, and
  the helper refuses one with the fix in the message. Assert the shape separately in the test, where a
  reader can see it.
- **Name as `area/thing`**, so `reference.toml` sorts into groups and a diff stays readable.
- **Prefer several specific numbers to one summary.** The per-subcommunity `subdiv` vectors catch a
  change that redistributes diversity between subcommunities while preserving the metacommunity
  figure; the `metadiv` scalar alone does not. Both are blessed for exactly that reason.
- **Bless both scales.** They aggregate with *different* power-mean orders — opposite ones for the
  relative-entropy measures — so a subcommunity vector can move while its metacommunity summary does
  not.
- **Keep ordinary assertions alongside the blessed ones.** A blessed number tells you *something
  changed*; a property tells you *the answer is still right*. Re-blessing silences the first and must
  never be able to silence the second. Both files assert numbers equivalence, the `β̄ · ρ̄ = 1`
  reciprocal, `α = ᾱ / w`, and symmetry where symmetry is required.
- **Nothing random.** Build trees explicitly with `createnode!` / `createbranch!` rather than
  `rand(Nonultrametric(…))`; a random topology re-blesses to noise every run. There is deliberately no
  seeding here — a fixed structure is readable in the diff, a seed is not.

⚠️ The function is `blessed`, not `canonical`, because **`BioSequences` exports `canonical`** (the
canonical orientation of a k-mer). A file that loads it — `test_types.jl` does — would get an
ambiguity between the two rather than either, reported as a bare `UndefVarError`. EcoSISTEM has no
such clash and calls it `canonical`; do not rename this one back to match.

## Notes

- A value with no blessed counterpart reports as **Broken**, not failed: adding a canonical test
  should not break the build before you have blessed it.
- `reference.toml` is **merged**, not replaced, on blessing. A partial run therefore cannot silently
  delete the blessed values of files it did not execute.
- The default tolerance is tight (`rtol = 1e-8`). Widen it only where a result genuinely is not
  reproducible to more digits, and say why at the call site.
