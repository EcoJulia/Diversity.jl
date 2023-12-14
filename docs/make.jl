# SPDX-License-Identifier: BSD-2-Clause

using Pkg
"Diversity" ∈ [p.name for p in values(Pkg.dependencies())] &&
    Pkg.rm("Diversity")
Pkg.develop(path = joinpath(@__DIR__, ".."))

using Documenter
using Diversity
using Phylo

makedocs(modules = [Diversity,
             Diversity.Ecology, Diversity.Jost,
             Diversity.Hill,
             Diversity.ShortNames, Diversity.API],
         sitename = "Diversity.jl")

deploydocs(repo = "github.com/EcoJulia/Diversity.jl.git",
           devbranch = "dev",
           push_preview = true)
