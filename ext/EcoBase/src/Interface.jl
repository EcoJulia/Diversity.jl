using Compat

# Functions - most have to be implemented with the concrete type
occurrences(asm::AbstractAssemblage) = error("function not defined for this type")
view(asm::AbstractAssemblage) = error("function not defined for this type")
places(asm::AbstractAssemblage) = error("function not defined for this type")
things(asm::AbstractAssemblage) = error("function not defined for this type")

nplaces(plc::AbstractPlaces) = error("function not defined for this type")
placenames(plc::AbstractPlaces) = error("function not defined for this type")

nthings(thg::AbstractThings) = error("function not defined for this type")
thingnames(thg::AbstractThings) = error("function not defined for this type")

nzrows(a::AbstractMatrix) = LinearIndices(Compat.sum(a .> 0, dims=2))[findall(Compat.sum(a .> 0, dims=2) .> 0)]
nzcols(a::AbstractMatrix) = LinearIndices(Compat.sum(a .> 0, dims=1))[findall(Compat.sum(a .> 0, dims=1) .> 0)]
nnz(a::AbstractArray) = Compat.sum(a .> 0)

occurring(asm::AbstractAssemblage) = nzrows(occurrences(asm))
occupied(asm::AbstractAssemblage) = nzcols(occurrences(asm))
if VERSION < v"0.7.0-"
    occupied(asm::AbstractAssemblage, idx) = collect(zip(findn(occurrences(asm)[idx, :])))
    occurring(asm::AbstractAssemblage, idx) = collect(zip(findn(occurrences(asm)[:, idx])))
else
    occupied(asm::AbstractAssemblage, idx) = findall(!iszero, occurrences(asm)[idx, :])
    occurring(asm::AbstractAssemblage, idx) = findall(!iszero, occurrences(asm)[:, idx])
end
noccurring(x) = length(occurring(x))
noccupied(x) = length(occupied(x))
noccurring(x, idx) = length(occurring(x, idx))
noccupied(x, idx) = length(occupied(x, idx))

thingoccurrences(asm::AbstractAssemblage, idx) = view(occurrences(asm), idx, :)
placeoccurrences(asm::AbstractAssemblage, idx) = view(occurrences(asm), :, idx) # make certain that the view implementation also takes thing or place names

richness(asm::AbstractAssemblage{Bool, T, P}) where {T, P} = vec(Compat.sum(occurrences(asm), dims=1))
richness(asm::AbstractAssemblage) = vec(mapslices(nnz, occurrences(asm), dims=1))

occupancy(asm::AbstractAssemblage{Bool, T, P}) where {T, P} = vec(Compat.sum(occurrences(asm), dims=2))
occupancy(asm::AbstractAssemblage) = vec(mapslices(nnz, occurrences(asm), dims=2))

records(asm::AbstractAssemblage) = nnz(occurrences(asm))

cooccurring(asm::AbstractAssemblage, inds...) = cooccurring(asm, [inds...])
function cooccurring(asm::AbstractAssemblage, inds::AbstractVector)
    sub = view(asm, species = inds)
    richness(sub) .== nthings(sub)
end

function createsummaryline(vec::AbstractVector{<:AbstractString})
    linefunc(vec) = mapreduce(x -> x * ", ", *, vec[1:(end-1)]) * vec[end]
    length(vec) == 1 && return vec[1]
    length(vec) < 6 && return linefunc(vec)
    linefunc(vec[1:3]) * "..." * linefunc(vec[(end-1):end])
end


function show(io::IO, asm::T) where T <: AbstractAssemblage
    tn = createsummaryline(thingnames(asm))
    pn = createsummaryline(placenames(asm))
    println(io,
    """$T with $(nthings(asm)) things in $(nplaces(asm)) places

    Thing names:
    $(tn)

    Place names:
    $(pn)
    """)
end

nplaces(asm::AbstractAssemblage, args...) = nplaces(places(asm), args...)
placenames(asm::AbstractAssemblage, args...) = placenames(places(asm), args...)
nthings(asm::AbstractAssemblage, args...) = nthings(things(asm), args...)
thingnames(asm::AbstractAssemblage, args...) = thingnames(things(asm), args...)

# TODO:
# accessing cache

xmin(grd::AbstractGrid) = error("function not defined for this type")
ymin(grd::AbstractGrid) = error("function not defined for this type")
xcellsize(grd::AbstractGrid) = error("function not defined for this type")
ycellsize(grd::AbstractGrid) = error("function not defined for this type")
xcells(grd::AbstractGrid) = error("function not defined for this type")
ycells(grd::AbstractGrid) = error("function not defined for this type")
cellsize(grd::AbstractGrid) = xcellsize(grd), ycellsize(grd)
cells(grd::AbstractGrid) = xcells(grd), ycells(grd)
xrange(grd::AbstractGrid) = xmin(grd):xcellsize(grd):xmax(grd) #includes intermediary points
yrange(grd::AbstractGrid) = ymin(grd):ycellsize(grd):ymax(grd)
xmax(grd::AbstractGrid) = xmin(grd) + xcellsize(grd) * (xcells(grd) - 1)
ymax(grd::AbstractGrid) = ymin(grd) + ycellsize(grd) * (ycells(grd) - 1)


indices(grd::AbstractGrid, idx) = error("function not defined for this type") #Implement this in SpatialEcology!
coordinates(grd::AbstractGrid) = error("function not defined for this type")
