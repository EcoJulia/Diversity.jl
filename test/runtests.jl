using Diversity
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

include("EffectiveNumbers.jl")
include("Metacommunity.jl")
include("DiversityMeasure.jl")
include("GeneralisedDiversities.jl")
include("Hill.jl")
include("Jost.jl")
include("Ecology.jl")
