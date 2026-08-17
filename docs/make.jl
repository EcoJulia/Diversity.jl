# SPDX-License-Identifier: BSD-2-Clause

using Pkg
Pkg.resolve()

using Documenter

# The modules, page order and site name are shared with `test/extras_docs.jl`, which builds the
# manual too so that a dangling cross-reference fails there rather than only here.
include("config.jl")

# Note: GR needs to be told there is no display, or plotting fails on a CI runner.
get!(ENV, "GKSwstype", "100")

makedocs(modules = DOCS_MODULES,
         sitename = DOCS_SITENAME,
         pages = DOCS_PAGES)

deploydocs(repo = "github.com/EcoJulia/Diversity.jl.git",
           devbranch = "dev",
           push_preview = true)
