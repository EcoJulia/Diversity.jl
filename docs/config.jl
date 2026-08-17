# SPDX-License-Identifier: BSD-2-Clause
#
# The parts of the Documenter configuration that `docs/make.jl` and `test/extras_docs.jl` have to
# agree on, kept in one place rather than written out twice. `extras_docs.jl` builds the manual to
# check its cross-references, and a second copy of this that drifted would check a different site
# from the one that gets published.
#
# `Phylo` is loaded here rather than in either caller because the extension has to be active for the
# phylogenetics page to build.

using Diversity
using Phylo

# Modules whose docstrings the manual pulls in through `@autodocs`, and against which `@ref` targets
# are resolved.
const DOCS_MODULES = [Diversity,
    Diversity.Ecology, Diversity.Jost,
    Diversity.Hill,
    Diversity.ShortNames, Diversity.API]

# The page order is set explicitly - otherwise Documenter sorts alphabetically and buries
# `framework.md`. A new page must be added here or it will not appear in the navigation.
const DOCS_PAGES = ["Introduction" => "index.md",
    "The framework" => "framework.md",
    "Building a metacommunity" => "metacommunities.md",
    "Coming from vegan" => "vegan.md",
    "Phylogenetic diversity" => "phylogenetics.md",
    "Genetic diversity" => "genetics.md",
    "Diversity.Ecology" => "ecology.md",
    "Diversity.Hill" => "hill.md",
    "Diversity.Jost" => "jost.md",
    "Diversity.API" => "api.md"]

const DOCS_SITENAME = "Diversity.jl"
