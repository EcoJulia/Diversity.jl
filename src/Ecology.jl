# SPDX-License-Identifier: BSD-2-Clause

using Diversity
using Diversity.ShortNames
using Diversity.API
using LinearAlgebra
using DataFrames
using EcoBase
using EcoBase: AbstractAssemblage

"""
    generalisedrichness(level::DiversityLevel, proportions::AbstractArray,
                        Z::AbstractMatrix)
    generalisedrichness(level::DiversityLevel, proportions::AbstractArray,
                        sim::AbstractTypes)

Calculates species richness (diversity at q = 0) of a series of
columns representing subcommunity counts, allowing a similarity matrix
for the types / species.

# Arguments:

- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

# Returns:

- diversity (at ecosystem level) or diversities (of subcommunities)
"""
function generalisedrichness(level::DiversityLevel,
                             proportions::AbstractArray,
                             Z::AbstractMatrix =
                             Matrix(1.0I, size(proportions, 1),
                                    size(proportions, 1)))
    return generalisedrichness(level, proportions, GeneralTypes(Z))
end

function generalisedrichness(level::DiversityLevel,
                             proportions::AbstractArray,
                             sim::AbstractTypes)
    if (level == subcommunityDiversity)
        dm = ᾱ
    elseif (level == metacommunityDiversity)
        dm = Gamma
    else
        error("Can't calculate richness for $level")
    end
    gr = level(dm(Metacommunity(proportions, sim)), 0)
    gr[!, :measure] .= "Richness"
    return gr
end

"""
    richness(proportions::AbstractMatrix)

Calculates species richness (diversity at q = 0) of a series of
columns representing independent subcommunity counts.

#### Arguments:
- `proportions`: population proportions

#### Returns:
- diversities of subcommunities
"""
function richness(proportions::AbstractVecOrMat)
    gr = generalisedrichness(subcommunityDiversity, proportions,
                             UniqueTypes(size(proportions, 1)))
    gr[!, :diversity] .= Int.(round.(gr[!, :diversity]))
    return gr
end

function richness(asm::AbstractAssemblage)
    hassimilarity(asm) &&
        error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return richness(occurrences(asm))
end

"""
    generalisedshannon(level::DiversityLevel, proportions::AbstractArray,
                       Z::AbstractMatrix)
    generalisedshannon(level::DiversityLevel, proportions::AbstractArray,
                       sim::AbstractTypes)

Calculates Shannon entropy (log of diversity at q = 1) of a series of
columns representing independent subcommunity counts, allowing a
similarity matrix for the types / species.

# Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

#### Returns:
- entropy (at metacommunity level) or entropies (of subcommunities)
"""
function generalisedshannon end

function generalisedshannon(level::DiversityLevel,
                            proportions::AbstractArray,
                            Z::AbstractMatrix = Matrix(1.0I,
                                                       size(proportions, 1),
                                                       size(proportions, 1)))
    return generalisedshannon(level, proportions, GeneralTypes(Z))
end

function generalisedshannon(level::DiversityLevel,
                            proportions::AbstractArray,
                            sim::AbstractTypes)
    if (level == subcommunityDiversity)
        dm = ᾱ
    elseif (level == metacommunityDiversity)
        dm = Gamma
    else
        error("Can't calculate richness for $level")
    end
    gs = level(dm(Metacommunity(proportions, sim)), 1)
    gs[!, :diversity] .= log.(gs[!, :diversity])
    gs[!, :measure] .= "Shannon"
    select!(gs, Not(:q))
    return gs
end

"""
    shannon(proportions::AbstractVecOrMat)

Calculates shannon entropy (log of diversity at q = 1) of a series of
columns representing independent subcommunity counts.

#### Arguments:
- `proportions`: population proportions

#### Returns:
- entropies of subcommunities
"""
function shannon(proportions::AbstractVecOrMat)
    return generalisedshannon(subcommunityDiversity,
                              proportions,
                              UniqueTypes(size(proportions,
                                               1)))
end

function shannon(asm::AbstractAssemblage)
    hassimilarity(asm) &&
        error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return shannon(occurrences(asm))
end

"""
    generalisedsimpson(level::DiversityLevel, proportions::AbstractArray,
                       Z::AbstractMatrix)
    generalisedsimpson(level::DiversityLevel, proportions::AbstractArray,
                       sim::AbstractTypes)

Calculates Simpson's index (1 / diversity at q = 2) of a series of
columns representing independent subcommunity counts, allowing a
similarity matrix for the types / species.

#### Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

#### Returns:
- concentration (at ecosystem level) or concentrations (of subcommunities)
"""
function generalisedsimpson end

function generalisedsimpson(level::DiversityLevel,
                            proportions::AbstractArray,
                            Z::AbstractMatrix = Matrix(1.0I,
                                                       size(proportions, 1),
                                                       size(proportions, 1)))
    return generalisedsimpson(level, proportions, GeneralTypes(Z))
end

function generalisedsimpson(level::DiversityLevel,
                            proportions::AbstractArray,
                            sim::AbstractTypes)
    if (level == subcommunityDiversity)
        dm = ᾱ
    elseif (level == metacommunityDiversity)
        dm = Gamma
    else
        error("Can't calculate richness for $level")
    end
    gs = level(dm(Metacommunity(proportions, sim)), 2)
    gs[!, :diversity] .= gs[!, :diversity] .^ -1
    gs[!, :measure] .= "Simpson"
    select!(gs, Not(:q))
    return gs
end

"""
    simpson(proportions::AbstractMatrix)

Calculates Simpson's index (1 / diversity at q = 2) of a series of
columns representing independent subcommunity counts.

#### Arguments:

- `proportions`: population proportions

#### Returns:

- concentrations of subcommunities
"""
function simpson(proportions::AbstractVecOrMat)
    return generalisedsimpson(subcommunityDiversity,
                              proportions,
                              UniqueTypes(size(proportions,
                                               1)))
end

function simpson(asm::AbstractAssemblage)
    hassimilarity(asm) &&
        error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return simpson(occurrences(asm))
end

"""
    generalisedjaccard(proportions::AbstractArray, qs, Z::AbstractMatrix)
    generalisedjaccard(proportions::AbstractArray, qs, sim::AbstractTypes)
    generalisedjaccard(meta::AbstractAssemblage, qs)

Calculates a generalisation of the Jaccard similarity of two columns
representing the counts of two subcommunities. This evaluates to raw
alpha / gamma - 1 for a series of orders, repesented as a vector of qs
(or a single number). It also includes an optional similarity matrix
for the species. This gives a measure of the distinctness of the
subcommunities, though we believe that beta and normalised beta have
better properties.

#### Arguments:

- `proportions`: population proportions
- `meta`: metacommunity / assemblage

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

#### Returns:

- Jaccard-related distinctivess measures
"""
function generalisedjaccard end

function generalisedjaccard(proportions::AbstractMatrix, Z::AbstractMatrix)
    return generalisedjaccard(proportions, GeneralTypes(Z))
end

function generalisedjaccard(proportions::AbstractMatrix, sim::AbstractTypes)
    return generalisedjaccard(Metacommunity(proportions, sim))
end

function generalisedjaccard(meta::AbstractAssemblage)
    countsubcommunities(meta) == 2 ||
        error("Can only calculate Jaccard index for 2 subcommunities")
    num = sum(minimum(getordinariness!(meta), dims = 2))
    denom = sum(maximum(getordinariness!(meta), dims = 2))
    jac = metadiv(Gamma(meta), 0)
    jac[!, :diversity] .= num / denom
    jac[!, :measure] .= "Jaccard"
    select!(jac, Not([:q]))
    return jac
end

"""
    jaccard(proportions::AbstractMatrix)
    jaccard(asm::AbstractAssemblage)

Calculates Jaccard similarity coefficient of two
columns representing independent subcommunity counts

#### Arguments:

- `proportions`: population proportions
- `asm`: assemblage / metacommunity

#### Returns:

- the Jaccard index
"""
jaccard(proportions::AbstractMatrix) = jaccard(Metacommunity(proportions))

function jaccard(asm::AbstractAssemblage)
    hassimilarity(asm) &&
        error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return generalisedjaccard(asm)
end

#=
"""
    generalisedpielou(level::DiversityLevel,
                      proportions::AbstractArray,
                      Z::AbstractMatrix)
    generalisedpielou(level::DiversityLevel,
                      proportions::AbstractArray[,
                      sim::AbstractTypes])
    generalisedpielou(level::DiversityLevel,
                      asm::AbstractAssemblage)

Calculates a generalisation of Pielou's evenness for columns
representing the counts or proportions of subcommunities. Values range from 
zero to one, with one representing complete evenness within the 
community. Since this is calculated as H / Hmax, this uses Shannon entropy
and q is effectively 1.

#### Arguments:
- `level`: DiversityLevel to calculate at (e.g. subcommunityDiversity)

- `proportions`: population proportions

- `Z`: similarity matrix or
- `sim`: instance of AbstractTypes

#### Returns:
- Pielou's evenness metric (at metacommunity level) or metrics (of subcommunities)
"""
=#
function generalisedpielou end
function generalisedpielou(level::DiversityLevel,
                           proportions::AbstractArray,
                           Z::AbstractMatrix)
    return generalisedpielou(level, Metacommunity(proportions, Z))
end

function generalisedpielou(level::DiversityLevel,
                           proportions::AbstractArray,
                           sim::AbstractTypes = UniqueTypes(size(proportions,
                                                                 1)))
    return generalisedpielou(level, Metacommunity(proportions, sim))
end

function generalisedpielou(level::DiversityLevel,
                           mc::AbstractAssemblage)
    hassimilarity(mc) &&
        error("Can't calculate Pielou evenness for $(typeof(gettypes(mc))) type as ill-defined maximum entropy")

    if (level == subcommunityDiversity)
        dm = ᾱ
        hmax = log.(level(dm(mc), 0).diversity)
    elseif (level == metacommunityDiversity)
        dm = Gamma
        mcab = Float64[x > 0 for x in getmetaabundance(mc)]
        mcab ./= sum(mcab)
        mc1 = Metacommunity(mcab, gettypes(mc))
        hmax = first(log.(level(dm(mc1), 0).diversity))
    else
        error("Can't calculate Pielou for $level")
    end
    gp = level(dm(mc), 1)
    gp[!, :diversity] .= log.(gp[!, :diversity]) ./ hmax
    gp[!, :measure] .= "Pielou"
    select!(gp, Not(:q))
    return gp
end

"""
    pielou(proportions::AbstractMatrix)
    pielou(asm::AbstractAssemblage)

Calculates Pielou's evenness of a series of
columns representing independent subcommunity counts.

#### Arguments:

- `proportions`: population proportions

#### Returns:

- evenness of subcommunities

#### Example:
```
communitymat = [10 20 30 20 0;
                10 0 50 80 10;
                60 10 90 0 0; 
                10 10 10 10 10;
                70 70 70 70 70;
                10 0 0 90 0];

pielou(communitymat)
```
"""
function pielou(proportions::AbstractVecOrMat)
    return generalisedpielou(subcommunityDiversity,
                             Metacommunity(proportions))
end

function pielou(asm::AbstractAssemblage)
    hassimilarity(asm) &&
        error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return generalisedpielou(subcommunityDiversity, asm)
end

"""
    gower(proportions::AbstractMatrix; countzeros::Bool = false, logscale::Bool = true)
    gower(asm::AbstractAssemblage; countzeros::Bool = false, logscale::Bool = true)

Calculates Gower's dissimilarity between exactly two subcommunities. Unlike the other measures here
this is a genuinely *pairwise* index, so it needs a metacommunity of two subcommunities and no more.

#### Arguments:

- `proportions`: population proportions; or
- `asm`: an AbstractAssemblage of exactly two subcommunities

#### Keyword arguments:

- `countzeros`: which of the two published conventions to use, and **not a free choice** - it
  selects the denominator. `true` divides by the total number of types, so types absent from both
  subcommunities still count; this is classic Gower (1971), and matches R vegan's `gower`. `false`
  divides only by the types actually present, matching vegan's `altGower`. Defaults to `false`.
- `logscale`: take `log10` of the abundances first, so that differences are proportional rather than
  absolute. Defaults to `false`.
- `normalise`: count each type as differing or not, rather than by how much it differs - a
  presence/absence reading. Defaults to whatever `countzeros` is.

#### Returns:

- Gower dissimilarity of the two subcommunities, as a single-row DataFrame

Note: R vegan 2.7 changed its `gower` to range-standardise columns first and drop tied columns from the
denominator, which returns `NA` for identical samples. This package keeps the classic reading, so
`gower(countzeros = true)` matches *old* vegan; see `test/run_rcall.jl`.
"""
function gower end

function gower(proportions::AbstractArray; countzeros::Bool = false,
               logscale::Bool = false, normalise::Bool = countzeros)
    return gower(Metacommunity(proportions), countzeros = countzeros,
                 logscale = logscale, normalise = normalise)
end

function gower(asm::AbstractAssemblage; countzeros::Bool = false,
               logscale::Bool = false, normalise::Bool = countzeros)
    countsubcommunities(asm) == 2 ||
        error("Can only calculate Gower distances for 2 subcommunities")

    g = meta_gamma(asm, 0)
    nz = countzeros ? nthings(asm) : noccurring(asm)
    occ = logscale ? [x == 0 ? 0 : log10(x) for x in getabundance(asm, true)] :
          getabundance(asm, true)
    diff = normalise ?
           sum(x > 0 ? 1 : 0 for x in abs.(occ[:, 1] .- occ[:, 2])) :
           sum(abs.(occ[:, 1] .- occ[:, 2]))
    g[!, :diversity] .= diff / nz
    g[!, :measure] .= "Gower"
    g[!, :countzeros] .= countzeros
    g[!, :logscale] .= logscale
    g[!, :normalise] .= normalise
    select!(g, Not(:q))
    return g
end

"""
    PAIRWISE_FORMS

The forms of [`pairwisedissimilarity`](@ref): `:jaccard`, `:sorensen`,
`:ochiai`, `:sokalsneath` and `:rao`.
"""
const PAIRWISE_FORMS = (:jaccard, :sorensen, :ochiai, :sokalsneath, :rao)

# The similarity of two subcommunities x and y in each form, from G = Pᵀ Z P
_pairwisesim(::Val{:jaccard}, G, gx, gy) = G / (gx + gy - G)
_pairwisesim(::Val{:sorensen}, G, gx, gy) = 2G / (gx + gy)
_pairwisesim(::Val{:ochiai}, G, gx, gy) = G / sqrt(gx * gy)
_pairwisesim(::Val{:sokalsneath}, G, gx, gy) = G / (2gx + 2gy - 3G)

"""
    pairwisedissimilarity(asm::AbstractAssemblage; forms = PAIRWISE_FORMS,
                          relative::Bool = false, squareroot::Bool = false)
    pairwisedissimilarity(proportions::AbstractMatrix[, Z::AbstractMatrix];
                          kwargs...)
    pairwisedissimilarity(proportions::AbstractMatrix, sim::AbstractTypes;
                          kwargs...)

Similarity-sensitive dissimilarity between every pair of subcommunities of
`asm`, for each of `forms`. Unlike [`jaccard`](@ref) and [`gower`](@ref) this is
not limited to two subcommunities, and it takes the similarity between types
into account, so it works with any types that have one (`GeneralTypes`,
`PhyloBranches`, genetic types, ...).

Every form is a function of G = Pᵀ Z P, the type similarity `Z` summed over
the pairs of individuals of two subcommunities (P holds the abundances of the
subcommunities as columns), which comes from the ordinariness Z P that the
metacommunity already caches. With g = diag(G), the dissimilarity between
subcommunities x and y is d = 1 - s for

  - `:jaccard`: s = G_xy / (g_x + g_y - G_xy)
  - `:sorensen`: s = 2 G_xy / (g_x + g_y)
  - `:ochiai`: s = G_xy / √(g_x g_y)
  - `:sokalsneath`: s = G_xy / (2g_x + 2g_y - 3G_xy)

(the similarity-based indices of Pavoine & Ricotta (2014), R adiv's
`dsimcom`), and, always on relative abundances, `:rao`: Rao's dissimilarity
coefficient (DISC) d = D_xy - ½(D_xx + D_yy), where D_xy = p_xᵀ Δ p_y is the
mean dissimilarity Δ = 1 - Z between an individual of x and one of y, which is
½(g_x + g_y) - G_xy.

With `Z` the identity (`UniqueTypes`) and presence/absence data the first four
are the classic Jaccard, Sørensen, Ochiai and Sokal-Sneath dissimilarities;
with relative abundances Sørensen is Morisita-Horn.

#### Arguments:

- `asm`: an assemblage or metacommunity, of two or more subcommunities; or
- `proportions`: the abundances of the subcommunities as columns, and
- `Z` or `sim`: the similarity of the types (default: none, `UniqueTypes`)

#### Keyword arguments:

- `forms`: which dissimilarities to calculate, from [`PAIRWISE_FORMS`](@ref)
- `relative`: use the relative abundances within each subcommunity rather than
  the abundances as given, so that each subcommunity has the same total weight.
  Defaults to `false`. `:rao` always uses relative abundances.
- `squareroot`: return √d rather than d. For positive semi-definite `Z` this
  makes each form Euclidean-embeddable (Pavoine & Ricotta 2014). Defaults to
  `false`.

#### Returns:

- a NamedTuple keyed by form of symmetric dissimilarity matrices over the
  subcommunities, in the order of `getsubcommunitynames(asm)`. Pairs involving
  an empty subcommunity are `NaN`.

A form that falls outside [0, 1] by rounding is clamped. A `Z` that is not
positive semi-definite can push a form further out; that is clamped too, with a
warning.

#### References:

- Pavoine, S. & Ricotta, C. (2014) Functional and phylogenetic similarity among
  communities. Methods in Ecology and Evolution 5, 666-675.
- Rao, C. R. (1982) Diversity and dissimilarity coefficients: a unified
  approach. Theoretical Population Biology 21, 24-43.
"""
function pairwisedissimilarity end

function pairwisedissimilarity(proportions::AbstractMatrix; kwargs...)
    return pairwisedissimilarity(Metacommunity(proportions); kwargs...)
end

function pairwisedissimilarity(proportions::AbstractMatrix, Z::AbstractMatrix;
                               kwargs...)
    return pairwisedissimilarity(proportions, GeneralTypes(Z); kwargs...)
end

function pairwisedissimilarity(proportions::AbstractMatrix, sim::AbstractTypes;
                               kwargs...)
    return pairwisedissimilarity(Metacommunity(proportions, sim); kwargs...)
end

function pairwisedissimilarity(asm::AbstractAssemblage;
                               forms = PAIRWISE_FORMS, relative::Bool = false,
                               squareroot::Bool = false)
    all(in(PAIRWISE_FORMS), forms) ||
        throw(ArgumentError("forms must be among $PAIRWISE_FORMS, got $forms"))
    n = countsubcommunities(asm)
    n ≥ 2 || error("Need at least two subcommunities, got $n")

    G = getabundance(asm, false)' * getordinariness!(asm)
    isapprox(G, G', rtol = sqrt(eps(Float64))) ||
        throw(ArgumentError("the similarity of " *
                            "$(typeof(gettypes(asm))) is not symmetric"))
    G = (G + G') / 2
    w = getweight(asm)
    Grel = G ./ (w * w')
    occupied = w .> 0

    ds = map(forms) do form
        Gf = form === :rao || relative ? Grel : G
        g = diag(Gf)
        D = fill(NaN, n, n)
        excess = 0.0
        for j in 1:n, i in 1:j
            (occupied[i] && occupied[j]) || continue
            if i == j
                D[i, i] = 0
                continue
            end
            d = form === :rao ? (g[i] + g[j]) / 2 - Gf[i, j] :
                1 - _pairwisesim(Val(form), Gf[i, j], g[i], g[j])
            excess = max(excess, -d, d - 1)
            d = clamp(d, 0, 1)
            D[i, j] = D[j, i] = squareroot ? sqrt(d) : d
        end
        excess > sqrt(eps(Float64)) &&
            @warn "$form: dissimilarity outside [0, 1] by up to $excess " *
                  "before clamping (is Z positive semi-definite?)"
        return D
    end
    return NamedTuple{Tuple(forms)}(ds)
end

"""
    generalisedfaith_pd(level::DiversityLevel, mc::AbstractMetacommunity)

Calculates Faith's phylogenetic diversity (PD) - the total length of the
branches of the tree spanned by the types present - at the level requested.

Unlike the other measures here, this one needs a tree, and is defined only for a
metacommunity built over the `PhyloBranches` supplied by the `Phylo` extension -
not for phylogenetic types in general, since the calculation relies on how that
type in particular maps leaf abundances onto branches. It is therefore available
only once `Phylo` is loaded.

It also takes no `q`: Faith's PD is the q = 0 case by definition, so there is no
profile to ask for, and PD depends on which types are present and not at all on
their abundances.

#### Arguments:
- `level`: the diversity level to calculate at - `subcommunityDiversity` for the
  PD of each subcommunity in isolation, `metacommunityDiversity` for the PD of
  the whole tree spanned by the metacommunity
- `mc`: a metacommunity built over a `PhyloBranches`

#### Returns:
- PD of the metacommunity, or PDs of the subcommunities
"""
function generalisedfaith_pd end

"""
    faith_pd(mc::AbstractMetacommunity)

Calculates Faith's phylogenetic diversity (PD) of each subcommunity of `mc` in
isolation - the total length of the branches of the tree spanned by the types
present in it. See [`generalisedfaith_pd`](@ref) for the metacommunity as a
whole, for why there is no `q` argument, and for why this needs a `PhyloBranches`
specifically.

#### Arguments:
- `mc`: a metacommunity built over a `PhyloBranches`

#### Returns:
- PDs of the subcommunities
"""
function faith_pd end
