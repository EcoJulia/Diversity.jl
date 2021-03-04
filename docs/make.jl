using Documenter
using Diversity

makedocs(modules = [Diversity,
                    Diversity.Ecology, Diversity.Jost, Diversity.Hill, Diversity.Phylogenetics, Diversity.API],
         sitename = "Diversity.jl")

deploydocs(deps = Deps.pip("pygments",
                           "mkdocs",
                           "mkdocs-material",
                           "python-markdown-math"))
