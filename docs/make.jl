# SPDX-License-Identifier: BSD-2-Clause

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
