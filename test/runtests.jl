using Diversity

tests = ["EffectiveNumbers", "Metacommunity", "DiversityMeasure",
         "GeneralisedDiversities", "Ecology", "Hill", "Jost"]

println("Running tests ...")

for t in tests
    fn = "test_$t.jl"
    println("* Testing $t.jl ...")
    include(fn)
end
