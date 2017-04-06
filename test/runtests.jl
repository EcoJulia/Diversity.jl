using Diversity

tests = ["EffectiveNumbers", "Metacommunity", "DiversityMeasure",
         "GeneralisedDiversities", "Ecology", "Hill", "Jost"]

println("Running tests ...")

for t in tests
    fn = "$t.jl"
    println("* $fn ...")
    include(fn)
end
