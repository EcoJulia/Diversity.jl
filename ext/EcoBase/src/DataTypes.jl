"""
    AbstractThings

Supertype for container of objects being observed, whether these are
species, sequences, tips of a phylogeny (which could be either), or
some other type of thing. This will contain the names of the things
being observed, and (optionally) metadata about them, such as a
phylogeny that connects them, taxonomic information, their sequences,
trait information, information on similarity between the different
things, etc.

"""
abstract type AbstractThings end

"""
    AbstractPlaces

Supertype for container of places where things are found (see
AbstractThings). This will contain names or a reference for the
places, and (optionally) metadata such as what kind of place these
are, e.g. geographic locations (see subtype AbstractLocations), where
additional metadata will includes location data, animal ids for where
the samples of things where retrieved from, etc.

"""
abstract type AbstractPlaces end

"""
    AbstractLocations <: AbstractPlaces

Subtype of AbstractPlaces that refers to locations with some
geographical component. This may be a series of arbitrarily arranged
points, a series of areas, or even grid of of regularly spaced
quadrats (see AbstractGrid).

"""
abstract type AbstractLocations <: AbstractPlaces end

"""
    AbstractGrid <: AbstractLocations <: AbstractPlaces

Subtype of AbstractLocations that refers to a grid of regularly
spaced, identically shaped, locations.
"""
abstract type AbstractGrid <: AbstractLocations end

"""
    AbstractAssemblage{D <: Real (e.g. Int, Float64, Bool),
                       T <: AbstractThings,
                       P <: AbstractPlaces}

An assemblage of things recorded as being present in one or more
places. These may, for instance, be species counts in quadrats over a
regular grid, relative abundance of viral sequences in a group of
individuals, or presence-absence of genera over multiple islands.

"""
abstract type AbstractAssemblage{D <: Real,
                                 T <: AbstractThings,
                                 P <: AbstractPlaces} end

# depend explicitly on AxisArrays (and/or sparse) in some implementation?
