module ValidateRCall

# Environment variable to avoid boring R package builds
mustCrossvalidate = haskey(ENV, "JULIA_MUST_CROSSVALIDATE") && ENV["JULIA_MUST_CROSSVALIDATE"] == "1"

# Only run R on unix or when R is installed because JULIA_MUST_CROSSVALIDATE is set to 1
global skipR = !mustCrossvalidate && !Sys.isunix()
try
    skipR && error("Skipping R testing...")
    using RCall
    global skipR = false
catch
    global skipR = true
    @warn "R or appropriate Phylo package not installed, skipping R cross-validation."
end

!skipR && include("run_rcall.jl")

end
