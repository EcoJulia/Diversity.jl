# Identify files in test/ that are testing matching files in src/
#  - src/Source.jl will be matched by test/test_Source.jl

# Need to install our version of PhyloTrees for the time-being
testphylo = true
if Pkg.installed("PhyloTrees") == nothing
    Pkg.clone("https://github.com/boydorr/PhyloTrees.jl.git")
elseif Pkg.installed("PhyloTrees") <= v"0.7.0"
    warn("Unable to run tests on phlyogenetic section")
    testphylo = false
end

filebase = map(file -> replace(file, r"(.*).jl", s"\1"),
                filter(file -> ismatch(r".*\.jl", file), readdir("../src")))
testbase = map(file -> replace(file, r"test_(.*).jl", s"\1"),
                filter(str -> ismatch(r"^test_.*\.jl$", str), readdir()))

info("Running tests for files:")
for t in testbase
    println("    = $t.jl")
end
println()

info("Running tests...")
for t in testbase
    if t != "Phylogenetics" || testphylo
        fn = "test_$t.jl"
        println("    * Testing $t.jl ...")
        include(fn)
        println()
    else
        fn = "test_$t.jl"
        println("    * Ignoring $t.jl as PhyloTrees not available...")
        println()
    end
end

# Identify tests with no matching file
superfluous = filter(f -> f ∉ filebase, testbase)
if length(superfluous) > 0
    info("Potentially superfluous tests:")
    for f in superfluous
        println("    + $f.jl")
    end
    println()
end

# Identify files with no matching test
missing = filter(f -> f ∉ testbase, filebase)
if length(missing) > 0
    info("Potentially missing tests:")
    for f in missing
        println("    - $f.jl")
    end
    println()
end

# Identify files that are cross-validating results against other packages
# test/pkg_Package.jl should validate results against the Package package

pkgbase = map(file -> replace(file, r"pkg_(.*).jl", s"\1"),
                   filter(str -> ismatch(r"^pkg_.*\.jl$", str), readdir()))

if length(pkgbase) > 0
    info("Running cross-validation against:")
    for p in pkgbase
        println("    = $p")
    end
    println()
    
    info("Running validation...")
    for p in pkgbase
        fn = "pkg_$p.jl"
        println("    * Validating against $p ...")
        include(fn)
        println()
    end
end
