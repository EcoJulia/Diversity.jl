# SPDX-License-Identifier: BSD-2-Clause

using Pkg

# Update Diversity folder packages 
Pkg.activate(".")
Pkg.update()

# Update examples folder packages
if isdir("examples")
    if isfile("examples/Project.toml")
        Pkg.activate("examples")
        Pkg.update()
        "Diversity" ∈ [p.name for p in values(Pkg.dependencies())] &&
            Pkg.rm("Diversity")
        Pkg.develop("Diversity")
    end
end

# Update docs folder packages
Pkg.activate("docs")
Pkg.update()
"Diversity" ∈ [p.name for p in values(Pkg.dependencies())] &&
    Pkg.rm("Diversity")
Pkg.develop("Diversity")

# Reformat files in package
using JuliaFormatter
using Diversity
format(Diversity)

# Carry out crosswalk for metadata
using ResearchSoftwareMetadata
ResearchSoftwareMetadata.crosswalk()
