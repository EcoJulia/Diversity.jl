using Documenter
using Diversity

makedocs(modules = [Diversity,
                    Diversity.Ecology, Diversity.Jost, Diversity.Hill, Diversity.ShortNames, Diversity.API],
         sitename = "Diversity.jl")

deploydocs(repo = "github.com/richardreeve/Diversity.jl.git",
           devbranch = "dev",
           deps = Deps.pip("pygments",
                           "mkdocs",
                           "mkdocs-material",
                           "python-markdown-math"))
