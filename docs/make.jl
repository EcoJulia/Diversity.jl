# SPDX-License-Identifier: BSD-2-Clause

using Pkg
"Diversity" ∈ [p.name for p in values(Pkg.dependencies())] &&
    Pkg.rm("Diversity")
Pkg.develop(path = joinpath(@__DIR__, ".."))

using Documenter
using Diversity
using Phylo

# ⚠️ The page order is set explicitly. Left to itself Documenter sorts alphabetically, which buries
# framework.md — the page that explains what every other page assumes — between ecology and genetics.
# The order here is: what the package is, what it measures, then the capabilities, the submodules of
# older measures, and finally the extension API.
makedocs(modules = [Diversity,
             Diversity.Ecology, Diversity.Jost,
             Diversity.Hill,
             Diversity.ShortNames, Diversity.API],
         sitename = "Diversity.jl",
         pages = ["Introduction" => "index.md",
             "The framework" => "framework.md",
             "Phylogenetic diversity" => "phylogenetics.md",
             "Genetic diversity" => "genetics.md",
             "Diversity.Ecology" => "ecology.md",
             "Diversity.Hill" => "hill.md",
             "Diversity.Jost" => "jost.md",
             "Diversity.API" => "api.md"])

deploydocs(repo = "github.com/EcoJulia/Diversity.jl.git",
           devbranch = "dev",
           push_preview = true)
