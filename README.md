# Diversity
[![Diversity](http://pkg.julialang.org/badges/Diversity_0.4.svg)](http://pkg.julialang.org/?pkg=Diversity&ver=0.4)
[![Diversity](http://pkg.julialang.org/badges/Diversity_0.5.svg)](http://pkg.julialang.org/?pkg=Diversity&ver=0.5)

[![Build status](https://travis-ci.org/richardreeve/Diversity.jl.svg?branch=master)](https://travis-ci.org/richardreeve/Diversity.jl?branch=master)
[![Build status (Windows)](https://ci.appveyor.com/api/projects/status/github/richardreeve/Diversity.jl?svg=true&branch=master)](https://ci.appveyor.com/project/richardreeve/diversity-jl/branch/master)
[![Coveralls status](https://img.shields.io/coveralls/richardreeve/Diversity.jl.svg)](https://coveralls.io/r/richardreeve/Diversity.jl?branch=master)
[![Codecov status](https://codecov.io/gh/richardreeve/Diversity.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/richardreeve/Diversity.jl)

Docs:
[![Stable Documentation](https://readthedocs.org/projects/diversityjl/badge/?version=stable)](http://diversityjl.readthedocs.org/en/stable/diversity/)
[![Latest Documentation](https://readthedocs.org/projects/diversityjl/badge/?version=latest)](http://diversityjl.readthedocs.org/en/latest/diversity/)

New Docs:
[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://richardreeve.github.io/Diversity.jl/stable)
[![Latest Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://richardreeve.github.io/Diversity.jl/latest)

**Diversity** is a [Julia](http://www.julialang.org) package that
provides functionality for measuring alpha, beta and gamma diversity
of supercommunities (e.g. ecosystems) and their constituent
subcommunities. It uses the diversity measures described in the arXiv
paper [arXiv:1404.6520 (q-bio.QM)](http://arxiv.org/abs/1404.6520),
*How to partition diversity*. It also provides a series of other
related and older diversity measures through sub-modules. Currently
these are all ecological diversity measures, but this will be expanded
through interfacing to BioJulia.

This package is still in alpha and under heavy development, and so we
do not guarantee its correctness, although we are aware of no problems
with it. Please
[raise an issue](https://github.com/richardreeve/Diversity.jl/issues)
if you find any problems.

Version 0.3, which will be the next release, will have significant
changes to the standard interface for calculating diversity and to the
output format to provide consistency with our R package
[rdiversity](https://github.com/boydorr/rdiversity). Older interfaces
will be deprecated. The new interface is currently available on master
if you want to try it out - it passes all of our tests, but currently
has a very clunky output format, which will be fixed next.

## Install

*Diversity* is in `METADATA` and can be installed via `Pkg.add("Diversity")`.

## Usage

#### Diversity

The main package provides basic numbers-equivalent diversity measures
(described in [Hill, 1973](http://www.jstor.org/stable/1934352)),
similarity-sensitive diversity measures (generalised from Hill, and
described in
[Leinster and Cobbold, 2012](http://www.esajournals.org/doi/abs/10.1890/10-2402.1)),
and related alpha, beta and gamma diversity measures at the level of
the supercommunity and its component subcommunities (generalised in
turn from Leinster and Cobbold, and described in
[Reeve et al, 2014](http://arxiv.org/abs/1404.6520)). The diversity
functions exist both with unicode names (e.g. ᾱ()), which are not
automatically exported as we feel they are too short and with matching
ascii names (e.g. NormalisedAlpha()), which are. We also provide a
general function for extract any diversity measure for a series of
subcommunity relative abundances. The full documentation can be found
[here](http://diversityjl.readthedocs.org/en/stable/diversity/).

Accessing the main functionality in the package is simple:

```julia_skip
using Diversity
...
diversities = superdiv(NormalisedAlpha(Supercommunity(proportions, Z)), [0, 1, 2, Inf])
diversity = subdiv(RawRho(Supercommunity(proportions, Z)), 2)
```

The package also provides sub-modules with other diversity measures:

#### Diversity.Ecology

Many existing ecological diversity measures can be derived from our
diversity measures, and so we provide them in the Diversity.Ecology
submodule along with generalised versions of them that relate to our
general measures of alpha, beta and gamma diversity at subcommunity
and supercommunity levels. The generalisations of species richness,
Shannon entropy and Simpson's index are the only standard measures we
are aware of whose subcommunity components sum directly to the
corresponding supercommunity measure (although note that Simpson's
index decreases for increased diversity, so small components are more
diverse). Documentation for these diversity measures can be found
[here](http://diversityjl.readthedocs.org/en/stable/ecology/).

#### Diversity.Hill

[Hill numbers](http://www.jstor.org/stable/1934352) are found in the
Diversity.Hill sub-module.
Documentation for these diversity measures can be found
[here](http://diversityjl.readthedocs.org/en/stable/hill/).

#### Diversity.Jost

Lou Jost's
[diversity](http://dx.doi.org/10.1111/j.2006.0030-1299.14714.x)
[measures](http://www.esajournals.org/doi/abs/10.1890/06-1736.1) are
found in the Diversity.Jost sub-module.
Documentation for these diversity measures is
[here](http://diversityjl.readthedocs.org/en/stable/jost/).

## Documentation

Documentation is generated by the Base documentation in Julia and
online via the
[Documenter](https://github.com/JuliaDocs/Documenter.jl) package.

### Usage

Accessing the documentation in Julia is easy in v0.4 onwards:

```julia
using Diversity

# Returns any documentation for the qDZ function and all qDZ methods
?qDZ

# Returns the specific documentation for that qD method call
?qD([0.1, 0.2, 0.7], 2)
```

The documentation is also available online.

### Stable branch

The online documentation for the current stable branch is
[here](http://diversityjl.readthedocs.org/en/stable/diversity/), and
API docs start
[here](http://diversityjl.readthedocs.org/en/stable/api/Diversity/).

### Master branch

The online documentation for the latest master (unreleased) branch is
[here](http://diversityjl.readthedocs.org/en/latest/diversity/), and
API docs start
[here](http://diversityjl.readthedocs.org/en/latest/api/Diversity/).
