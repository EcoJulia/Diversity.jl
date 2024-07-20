# SPDX-License-Identifier: BSD-2-Clause

using Pkg
Pkg.update()
"Diversity" ∈ [p.name for p in values(Pkg.dependencies())] &&
    Pkg.rm("Diversity")
Pkg.develop("Diversity")

using Documenter
using Diversity
using Phylo

makedocs(modules = [Diversity,
             Diversity.Ecology, Diversity.Jost,
             Diversity.Hill,
             Diversity.ShortNames, Diversity.API],
         sitename = "Diversity.jl")

deploydocs(repo = "github.com/EcoJulia/Diversity.jl.git",
           push_preview = true,
           devbranch = "dev",
           devurl = "dev")
