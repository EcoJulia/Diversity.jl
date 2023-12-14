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
