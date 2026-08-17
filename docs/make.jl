# SPDX-License-Identifier: BSD-2-Clause

using Pkg
"Diversity" ∈ [p.name for p in values(Pkg.dependencies())] &&
    Pkg.rm("Diversity")
Pkg.develop(path = joinpath(@__DIR__, ".."))

using Documenter
using Diversity
using Phylo

# Note: GR needs to be told there is no display, or plotting fails on a CI runner.
get!(ENV, "GKSwstype", "100")

# The page order is set explicitly - otherwise Documenter sorts alphabetically.
makedocs(modules = [Diversity,
             Diversity.Ecology, Diversity.Jost,
             Diversity.Hill,
             Diversity.ShortNames, Diversity.API],
         sitename = "Diversity.jl",
         pages = ["Introduction" => "index.md",
             "The framework" => "framework.md",
             "Building a metacommunity" => "metacommunities.md",
             "Coming from vegan" => "vegan.md",
             "Phylogenetic diversity" => "phylogenetics.md",
             "Genetic diversity" => "genetics.md",
             "Diversity.Ecology" => "ecology.md",
             "Diversity.Hill" => "hill.md",
             "Diversity.Jost" => "jost.md",
             "Diversity.API" => "api.md"])

deploydocs(repo = "github.com/EcoJulia/Diversity.jl.git",
           devbranch = "dev",
           push_preview = true)
