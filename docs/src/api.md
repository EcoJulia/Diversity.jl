# Diversity.API

The **Diversity.API** submodule provides the API that must be extended
for new `AbstractTypes`, `AbstractPartition` and
`AbstractMetacommunity` subtypes. It is what lets a package that knows nothing
about this one — [EcoSISTEM](https://github.com/EcoJulia/EcoSISTEM.jl), for
instance — have its own types measured by every diversity measure here.

## Usage

Extending the system means adding methods to the underscore-prefixed functions in
`Diversity.API`. You add methods; you never replace them.

A new `AbstractTypes` subtype needs only two: `_gettypenames` and
`_calcsimilarity`. Everything else has a working default — `_counttypes`, for
example, falls back to the length of the type name vector. Here is a complete
one, in which every distinct pair of types has the same similarity:

```@repl api
using Diversity
using Diversity.API

struct UniformSimilarity <: Diversity.API.AbstractTypes
    names::Vector{String}
    similarity::Float64
end

Diversity.API._gettypenames(us::UniformSimilarity, ::Bool) = us.names

function Diversity.API._calcsimilarity(us::UniformSimilarity, ::Real)
    n = length(us.names)
    Z = fill(us.similarity, n, n)
    for i in 1:n
        Z[i, i] = 1.0
    end
    return Z
end
```

That is enough to use it everywhere a built-in type would go:

```@repl api
types = UniformSimilarity(["a", "b", "c"], 0.5)
counttypes(types)
meta = Metacommunity([0.5, 0.3, 0.2], types)
meta_gamma(meta, 0)
```

The similarity parameter shows what the measures are doing. At `0.0` no type
resembles any other, so the metacommunity holds three types' worth of diversity;
at `1.0` every type is interchangeable and it holds one:

```@repl api
meta_gamma(Metacommunity([0.5, 0.3, 0.2],
                         UniformSimilarity(["a", "b", "c"], 0.0)), 0)
meta_gamma(Metacommunity([0.5, 0.3, 0.2],
                         UniformSimilarity(["a", "b", "c"], 1.0)), 0)
```

## The contract

| for a new... | must implement | may implement — and what you get if you do not |
|---|---|---|
| `AbstractTypes` | `_gettypenames`, `_calcsimilarity` (unless `_hassimilarity` is `false`) | `_counttypes` (counts the type names), `_calcabundance` (the abundances unchanged, with scale `1`), `_calcordinariness` (`_calcsimilarity(t, scale) * abundances`), `_getdiversityname` (`"unknown"`), `_addedoutputcols` (no extra columns), `_getaddedoutput` (`nothing`), `floattypes` (every `AbstractFloat`), `_hassimilarity` (`true`), `_subsettypes` (materialises the similarity into a `GeneralTypes`) |
| `AbstractPartition` | `_getsubcommunitynames` | `_countsubcommunities` (counts the subcommunity names), `_subsetpartition` (builds a `Subcommunities` from the names kept) |
| `AbstractMetacommunity` | `_gettypes`, `_getpartition`, `_getabundance` | `_getmetaabundance` (abundances summed across subcommunities), `_getweight` (abundances summed across types), `_getordinariness!` (`_calcordinariness` of the types, abundances and scale), `_getmetaordinariness!` (ordinariness summed across subcommunities), `_getscale` (`1`) |

Two of the optional ones are worth knowing about even if you do not implement
them. `_calcabundance` returns both the processed abundances and a **scale**,
which is then passed to `_calcsimilarity` — that is what lets a phylogeny measure
diversity over branches rather than over species. And the `raw::Bool` argument
carried through the API distinguishes the types the user supplied from the types
diversity is actually computed over, which differ for exactly that reason.

Leave a required method out and you get an error naming it. That is worth saying
because these abstract types are subtypes of `EcoBase`'s, so the fallbacks that
let a plain `EcoBase` assemblage be measured directly would otherwise call back
into the method you had not written, and the symptom would be a stack overflow
rather than a missing method.

`_calcsimilarity` is the one required method with a way out, and it is the pair
to `_hassimilarity` above. Declare `_hassimilarity(::YourTypes) = false` and you
need not write it: you get an identity matrix, every type like itself and nothing
else, which is what `UniqueTypes` means. Claim similarity and omit the matrix and
you get the error instead — otherwise that same identity matrix would be used for
you, and the diversities would come back quietly wrong rather than not at all.

```@contents
```

```@autodocs
Modules = [Diversity.API]
Private = false
```

```@index
```
