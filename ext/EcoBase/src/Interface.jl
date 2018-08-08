# Functions - most have to be implemented with the concrete type
occurrences(asm::AbstractAssemblage) = error("function not defined for this type")
view(asm::AbstractAssemblage) = error("function not defined for this type")
places(asm::AbstractAssemblage) = error("function not defined for this type")
things(asm::AbstractAssemblage) = error("function not defined for this type")

nplaces(plc::AbstractPlaces) = error("function not defined for this type")
placenames(plc::AbstractPlaces) = error("function not defined for this type")

nthings(thg::AbstractThings) = error("function not defined for this type")
thingnames(thg::AbstractThings) = error("function not defined for this type")

nzrows(a::AbstractMatrix) = find(sum(a .> 0, 2) .> 0)
nzcols(a::AbstractMatrix) = find(sum(a .> 0, 1) .> 0)
nnz(a::AbstractArray) = sum(a .> 0)

occurring(asm::AbstractAssemblage) = nzrows(occurrences(asm))
occupied(asm::AbstractAssemblage) = nzcols(occurrences(asm))
occupied(asm::AbstractAssemblage, idx) = findn(occurrences(asm)[idx, :])
occurring(asm::AbstractAssemblage, idx) = findn(occurrences(asm)[:, idx])

noccurring(x) = length(occurring(x))
noccupied(x) = length(occupied(x))
noccurring(x, idx) = length(occurring(x, idx))
noccupied(x, idx) = length(occupied(x, idx))

thingoccurrences(asm::AbstractAssemblage, idx) = view(occurrences(asm), idx, :)
placeoccurrences(asm::AbstractAssemblage, idx) = view(occurrences(asm), :, idx) # make certain that the view implementation also takes thing or place names

richness(asm::AbstractAssemblage{Bool, T, P}) where {T, P} = vec(sum(occurrences(asm), 1))
richness(asm::AbstractAssemblage) = vec(mapslices(nnz, occurrences(asm), 1))

occupancy(asm::AbstractAssemblage{Bool, T, P}) where {T, P} = vec(sum(occurrences(asm), 2))
occupancy(asm::AbstractAssemblage) = vec(mapslices(nnz, occurrences(asm), 2))

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


macro forward_func(ex, fs)
    T, field = ex.args[1], ex.args[2].args[1]
    T = esc(T)
    fs = Meta.isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
    :($([:($f(x::$T, args...) = (Base.@_inline_meta; $f($(field)(x), args...)))
        for f in fs]...);
    nothing)
end

@forward_func AbstractAssemblage.places nplaces, placenames
@forward_func AbstractAssemblage.things nthings, thingnames

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
