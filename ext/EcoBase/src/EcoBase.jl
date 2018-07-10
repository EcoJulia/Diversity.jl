module EcoBase

import Base: show
import RecipesBase

include("DataTypes.jl")

include("Interface.jl")
export nthings, nplaces, occupancy, richness, records, placenames, thingnames
export occurring, noccurring, occupied, noccupied, occurrences
export indices, coordinates, xcells, ycells, cells, xmin, xmax, ymin, ymax
export xrange, yrange, xcellsize, ycellsize, cellsize

include("PlotRecipes.jl")
end # module
