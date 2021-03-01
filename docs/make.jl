using Documenter
using Diversity

makedocs(modules = [Diversity,
                    Diversity.Ecology, Diversity.Jost, Diversity.Hill],
         clean   = false,
         sitename = "Diversity.jl")

deploydocs(deps = Deps.pip("pygments",
                           "mkdocs",
                           "mkdocs-material",
                           "python-markdown-math"),
           repo = "github.com/richardreeve/Diversity.jl.git")
