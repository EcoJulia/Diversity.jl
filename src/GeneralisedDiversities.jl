# SPDX-License-Identifier: BSD-2-Clause

using EcoBase: AbstractAssemblage
using Diversity.ShortNames

# The fourteen wrappers below each name one measure at one scale, so that a caller who knows what
# they want to measure need not also know which DiversityMeasure and which of subdiv/metadiv
# produces it. Each is a thin call through to that pair; see the framework page of the documentation
# for what the measures mean.

"""
    norm_sub_alpha(meta::AbstractAssemblage, qs)

Calculates the similarity-sensitive diversity of each subcommunity in isolation — the diversity
each subcommunity would have if it were the whole of the community.

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per subcommunity per order.
"""
function norm_sub_alpha(meta::AbstractAssemblage, qs)
    return subdiv(ᾱ(meta), qs)
end

"""
    raw_sub_alpha(meta::AbstractAssemblage, qs)

Calculates the per-subcommunity estimate of naive-community metacommunity diversity — the diversity
of the metacommunity that this subcommunity alone would imply, if no type were shared with any
other subcommunity. It is [`norm_sub_alpha`](@ref) per individual, rescaled by the size of the
subcommunity.

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per subcommunity per order.
"""
function raw_sub_alpha(meta::AbstractAssemblage, qs)
    return subdiv(α(meta), qs)
end

"""
    norm_sub_beta(meta::AbstractAssemblage, qs)

Calculates the per-subcommunity estimate of the effective number of distinct subcommunities. It is
high when a subcommunity is both distinctive and small, and is the reciprocal of
[`norm_sub_rho`](@ref).

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per subcommunity per order.
"""
function norm_sub_beta(meta::AbstractAssemblage, qs)
    return subdiv(β̄(meta), qs)
end

"""
    raw_sub_beta(meta::AbstractAssemblage, qs)

Calculates the distinctiveness of individual subcommunities — how much of each subcommunity is
unlike the rest of the metacommunity, whether through types found nowhere else or through low
similarity to the types that are. It is the reciprocal of [`raw_sub_rho`](@ref).

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per subcommunity per order.
"""
function raw_sub_beta(meta::AbstractAssemblage, qs)
    return subdiv(β(meta), qs)
end

"""
    norm_sub_rho(meta::AbstractAssemblage, qs)

Calculates the representativeness of individual subcommunities — how typical each subcommunity is
of the metacommunity as a whole. A subcommunity holding a fixed fraction of equally abundant types
has representativeness equal to that fraction.

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per subcommunity per order.
"""
function norm_sub_rho(meta::AbstractAssemblage, qs)
    return subdiv(ρ̄(meta), qs)
end

"""
    raw_sub_rho(meta::AbstractAssemblage, qs)

Calculates the redundancy of individual subcommunities — the extent to which the diversity of the
metacommunity would survive the loss of each subcommunity. It takes its minimum of 1 when nothing
resembling the subcommunity remains elsewhere.

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per subcommunity per order.
"""
function raw_sub_rho(meta::AbstractAssemblage, qs)
    return subdiv(ρ(meta), qs)
end

"""
    sub_gamma(meta::AbstractAssemblage, qs)

Calculates the contribution per individual in a subcommunity toward metacommunity diversity. It
combines the subcommunity's own diversity with the rarity of its types in the metacommunity, so a
subcommunity of a few very rare types contributes heavily however dull it looks in isolation.

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per subcommunity per order.
"""
function sub_gamma(meta::AbstractAssemblage, qs)
    return subdiv(Γ(meta), qs)
end

"""
    norm_meta_alpha(meta::AbstractAssemblage, qs)

Calculates the average similarity-sensitive diversity of the subcommunities, each taken in
isolation. It is invariant under shattering — subdividing a subcommunity into identical parts does
not change it.

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per order.
"""
function norm_meta_alpha(meta::AbstractAssemblage, qs)
    return metadiv(ᾱ(meta), qs)
end

"""
    raw_meta_alpha(meta::AbstractAssemblage, qs)

Calculates naive-community metacommunity diversity — the diversity the metacommunity would have if
its subcommunities shared no types and no similarity. It is an upper bound on the true
metacommunity diversity [`meta_gamma`](@ref).

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per order.
"""
function raw_meta_alpha(meta::AbstractAssemblage, qs)
    return metadiv(α(meta), qs)
end

"""
    norm_meta_beta(meta::AbstractAssemblage, qs)

Calculates the effective number of distinct subcommunities. It reaches its maximum, the number of
subcommunities, when they are completely distinct and of equal size, and is invariant under
shattering.

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per order.
"""
function norm_meta_beta(meta::AbstractAssemblage, qs)
    return metadiv(β̄(meta), qs)
end

"""
    raw_meta_beta(meta::AbstractAssemblage, qs)

Calculates the average distinctiveness of the subcommunities. It can be read as a kind of turnover
between each subcommunity and the rest of the metacommunity.

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per order.
"""
function raw_meta_beta(meta::AbstractAssemblage, qs)
    return metadiv(β(meta), qs)
end

"""
    norm_meta_rho(meta::AbstractAssemblage, qs)

Calculates the average representativeness of the subcommunities.

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per order.
"""
function norm_meta_rho(meta::AbstractAssemblage, qs)
    return metadiv(ρ̄(meta), qs)
end

"""
    raw_meta_rho(meta::AbstractAssemblage, qs)

Calculates the average redundancy of the subcommunities. It takes its minimum of 1 when the
subcommunities have nothing in common, and rises towards the number of subcommunities as they
become more alike.

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per order.
"""
function raw_meta_rho(meta::AbstractAssemblage, qs)
    return metadiv(ρ(meta), qs)
end

"""
    meta_gamma(meta::AbstractAssemblage, qs)

Calculates metacommunity similarity-sensitive diversity — the diversity of the metacommunity taken
as a whole, ignoring how it is divided. It is the average of the subcommunity contributions
[`sub_gamma`](@ref).

# Arguments:

- `meta`: a Metacommunity
- `qs`: a single order or a vector of orders

# Returns:

- A DataFrame of diversities, one row per order.
"""
function meta_gamma(meta::AbstractAssemblage, qs)
    return metadiv(Γ(meta), qs)
end

"""
### Calculates subcommunity and metacommunity diversities

Calculates any diversity of a Metacommunity for a series of orders,
repesented as one or a vector of qs.

#### Arguments:
- `dls`: an iterable collection of DiversityLevels
- `dms`: an iterable collection of DiversityMeasures
- `meta`: a Metacommunity
- `qs`: single number or vector of values of parameter q

#### Returns:

A vector containing all of the diversity levels of all of the requested diversities.
"""
function diversity(dls, dms, meta::AbstractAssemblage, qs)
    return mapreduce(measure -> mapreduce(dl -> dl(measure, qs), append!, dls),
                     append!, map(dm -> dm(meta), dms))
end
