using Documenter
using Diversity

makedocs(modules = [Diversity, Diversity.Phylogenetics,
                    Diversity.Ecology, Diversity.Jost, Diversity.Hill],
         clean   = false)

deploydocs(deps = Deps.pip("pygments",
                           "mkdocs",
                           "mkdocs-material",
                           "python-markdown-math"),
           repo = "github.com/richardreeve/Diversity.jl.git",
           julia="0.6",
           osname="linux")
