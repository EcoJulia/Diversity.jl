## Diversity

**Diversity** is a [Julia](http://www.julialang.org) package that
provides functionality for measuring alpha, beta and gamma diversity
of metacommunities (e.g. ecosystems) and their constituent
subcommunities. It uses the diversity measures described in the arXiv
paper [arXiv:1404.6520 (q-bio.QM)][paper-url],
*How to partition diversity*. It also provides a series of other
related and older diversity measures through sub-modules. Currently
these are all ecological diversity measures, but this will be
expanded, possibly through interfacing to BioJulia.

This package is in beta now, but is cross-validated against our R
package [boydorr/rdiversity][rdiversity-url], which is developed
independently, so please [raise an issue][issues-url] if you find any
problems. We now use a DataFrame as the common output format for all
of the diversity calculations to provide consistency with our R
package [rdiversity][rdiversity-url]. The code is certainly not
optimised for speed at the moment due to the substantial changes that
have happened to it under the hood, and the Phylogenetics submodule is
also new, and may need further improvements.

Version 0.4, which has been recently released, has significant
breaking changes to the underlying code, whuch mean it is no longer
compatible with Julia v0.5. It is periodically working with Julia
nightly and I aim to keep it that way (though other packages need to
be updated too). Older interfaces from v0.2 have been removed in v0.4.

## Install

*Diversity* is in `METADATA` and can be installed via `Pkg.add("Diversity")`.

## Usage

### Diversity

The main package provides basic numbers-equivalent diversity measures
(described in [Hill, 1973](http://www.jstor.org/stable/1934352)),
similarity-sensitive diversity measures (generalised from Hill, and
described in
[Leinster and Cobbold, 2012][leinster-cobbold-url]),
and related alpha, beta and gamma diversity measures at the level of
the metacommunity and its component subcommunities (generalised in
turn from Leinster and Cobbold, and described in
[arXiv:1404.6520 (q-bio.QM)][paper-url]). The diversity functions
exist both with unicode names (e.g. ᾱ()), which are not automatically
exported as we feel they are too short and with matching ascii names
(e.g. NormalisedAlpha()), which are. We also provide a general
function for extract any diversity measure for a series of
subcommunity relative abundances.

#### Getting started

Before calculating diversity a `Metacommunity` object must be created. This object contains all the information needed to calculate diversity.

```julia
# Load the package into Julia
using Diversity

# Example population
pop = [1 1 0; 2 0 0; 3 1 4]
pop = pop / sum(pop)

# Create Metacommunity object
meta = Metacommunity(pop)
```

#### Calculating diversity
First we need to calculate the low-level diversity component seperately, by passing a `metacommunity` object to the appropriate function; `RawAlpha()`, `NormalisedAlpha()`, `RawBeta()`, `NormalisedBeta()`, `RawRho()`, `NormalisedRho()`, or `Gamma()`.

```julia
# First, calculate the normalised alpha component
component = NormalisedAlpha(meta)
```

Afterwhich, `subdiv()` or `metadiv()` are used to calculate subcommunity or metacommunity diversity, respectively (since both subcommunity and metacommunity diversity measures are transformations of the same low-level components, this is computationally more efficient).

```julia
# Then, calculate species richness of the subcommunities
subdiv(component, 0)

# or the average (alpha) species richness across the whole population
metadiv(component, 0)

# We can also generate a diversity profile by calculating multiple q-values simultaneously
df = subdiv(component, 0:30)
```

In some instances, it may be useful to calculate **all** subcommunity (or metacommunity) measures. In which case, a `Metacommunity` object may be passed directly to `subdiv()` or `metadiv()`:

```julia
# To calculate all subcommunity diversity measures
subdiv(meta, 0:2)

# To calculate all metacommunity diversity measures
metadiv(meta, 0:2)
```

Alternatively, if computational efficiency is not an issue, a single measure of diversity may be calculated directly by calling a wrapper function:
```julia
norm_sub_alpha(meta, 0:2)
```
A complete list of these functions is shown below:

* `raw_sub_alpha()` : per-subcommunity estimate of naive-community metacommunity diversity
* `norm_sub_alpha()` : similarity-sensitive diversity of each subcommunity in isolation
* `raw_sub_rho()` : redundancy of individual subcommunities
* `norm_sub_rho()` : representativeness of individual subcommunities
* `raw_sub_beta()` : distinctiveness of individual subcommunities
* `norm_sub_beta()` : per-subcommunity estimate of effective number of distinct subcommunities
* `sub_gamma()` : contribution per individual in a subcommunity toward metacommunity diversity
* `raw_meta_alpha()` : naive-community metacommunity diversity
* `norm_meta_alpha()` : average similarity-sensitive diversity of subcommunities
* `raw_meta_rho()` : average redundancy of subcommunities
* `norm_meta_rho()` : average representativeness of subcommunities
* `raw_meta_beta()` : average distinctiveness of subcommunities
* `norm_meta_beta()` : effective number of distinct subcommunities
* `meta_gamma()` : metacommunity similarity-sensitive diversity

#### Diversity.Phylogenetics

Phylogenetic diversity (described [here][leinster-cobbold-url]) is
included in the Diversity.Phylogenetics submodule. Documentation for
these diversity measures can be found
[here](http://richardreeve.github.io/Diversity.jl/latest/phylogenetics/).
The phylogenetics code relies on the [Phylo][phylo-url] package to
generate trees to incorporate into the diversity code:

```julia-repl
julia> using Diversity, Phylo, Diversity.Phylogenetics
Creating Diversity to Phylo interface...

julia> communities = [4 1; 3 2; 1 0; 0 1] / 12;

julia> nt = rand(Nonultrametric(4))
NamedTree phylogenetic tree with 7 nodes and 6 branches
Leaf names:
String["tip 1", "tip 2", "tip 3", "tip 4"]

julia> metaphylo = Metacommunity(communities, PhyloTypes(nt));

julia> raw_meta_rho(metaphylo, [1, 2])
2×8 DataFrames.DataFrame
│ Row │ div_type     │ measure │ q │ type_level │ type_name │ partition_level │ partition_name │ diversity │
├─────┼──────────────┼─────────┼───┼────────────┼───────────┼─────────────────┼────────────────┼───────────┤
│ 1   │ Phylogenetic │ RawRho  │ 1 │ types      │           │ metacommunity   │       │ 1.80666   │
│ 2   │ Phylogenetic │ RawRho  │ 2 │ types      │           │ metacommunity   │       │ 1.7056    │
```

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://richardreeve.github.io/Diversity.jl/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://richardreeve.github.io/Diversity.jl/stable

[travis-img]: https://travis-ci.org/richardreeve/Diversity.jl.svg?branch=master
[travis-url]: https://travis-ci.org/richardreeve/Diversity.jl?branch=master

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/richardreeve/Diversity.jl?svg=true&branch=master
[appveyor-url]: https://ci.appveyor.com/project/richardreeve/diversity-jl/branch/master

[coveralls-img]: https://img.shields.io/coveralls/richardreeve/Diversity.jl.svg
[coveralls-url]: https://coveralls.io/r/richardreeve/Diversity.jl?branch=master

[codecov-img]: https://codecov.io/gh/richardreeve/Diversity.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/richardreeve/Diversity.jl

[issues-url]: https://github.com/richardreeve/Diversity.jl/issues

[pkg-0.5-img]: http://pkg.julialang.org/badges/Diversity_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/?pkg=Diversity&ver=0.5
[pkg-0.6-img]: http://pkg.julialang.org/badges/Diversity_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=Diversity&ver=0.6

[paper-url]: http://arxiv.org/abs/1404.6520

[rdiversity-url]: https://github.com/boydorr/rdiversity

[leinster-cobbold-url]: http://www.esajournals.org/doi/abs/10.1890/10-2402.1

[phylo-url]: https://github.com/richardreeve/Phylo.jl
