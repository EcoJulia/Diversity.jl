using Diversity

"""
### Calculates Jost's alpha diversity

Calculates Jost's alpha diversity of a series of columns representing
independent community counts, for a series of orders, repesented as a
vector of qs. This is just the naive-community ecosystem diversity
divided by the naive-community beta diversity.

### Arguments:
- `proportions` relative proportions of different individuals / species
                in population (vector, or matrix where columns are
                for individual sub-communities)

- `qs` single number or vector of orders of diversity measurement

### Returns:
- array of diversities, first dimension representing sub-communities, and
  last representing values of q
"""
function jostalpha{S <: AbstractFloat}(proportions::Matrix{S}, qs)
    metacommunityDiversity(RawAlpha(Metacommunity(proportions)), qs) ./
    qD(reshape(mapslices(sum, proportions, (1,)), size(proportions)[2]), qs)
end

"""
### Calculates Jost's beta diversity

Calculates Jost's beta diversity of a series of columns representing
independent community counts, for a series of orders, repesented as a
vector of qs. This is just the naive gamma diversity divided by
Jost's alpha diversity

### Arguments:
- `proportions` relative proportions of different individuals / species
                in population (vector, or matrix where columns are
                for individual sub-communities)

- `qs` single number or vector of orders of diversity measurement

### Returns:
- array of diversities, first dimension representing sub-communities, and
  last representing values of q
"""
function jostbeta{S <: AbstractFloat}(proportions::Matrix{S}, qs)
    metacommunityDiversity(Gamma(Metacommunity(proportions)), qs) ./ jostalpha(proportions, qs)
end
