module TestPhylogenetics
using Base.Test
if !isdefined(Base.Test, Symbol("@test_warn"))
    # Ignore @test_warn unless it's there...
    macro test_warn(str, test)
    end
    macro test_nowarn(test)
    end
end

using Diversity
using DataFrames
using PhyloTrees
info("Testing phylogenetics")

end
