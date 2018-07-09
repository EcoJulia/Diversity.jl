module TestDiversity
using Compat.Test

using Diversity

@testset "Deprecations" begin
    # Only run deprecation checks on macs
    if is_apple()
        # No current deprecations?
    end
end

end
