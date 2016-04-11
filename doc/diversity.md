## Diversity

The main **Diversity** package provides basic numbers-equivalent
diversity measures (described in
[Hill, 1973](http://www.jstor.org/stable/1934352)),
similarity-sensitive diversity measures (generalised from Hill, and
described in
[Leinster and Cobbold, 2012](http://www.esajournals.org/doi/abs/10.1890/10-2402.1)),
and related alpha, beta and gamma diversity measures at the level of
the supercommunity and its component subcommunities (generalised in
turn from Leinster and Cobbold, and described in
[Reeve et al, 2014](http://arxiv.org/abs/1404.6520)). The diversity
functions exist both with unicode names (e.g. ```ᾱ()```), which are
not automatically exported (as we feel they are too short) and with
matching longer ASCII names (e.g. ```NormalisedAlpha()```), which are.
We also provide functions to calculate appropriate
```subcommunityDiversity()``` and ```supercommunityDiversity()```
values for each measure, a general ```diversity()``` function for
extract any diversity measure at a series of scales.

Accessing the main functionality in the package is simple:

```julia_skip
using Diversity
...
diversities = supercommunityDiversity(NormalisedAlpha(Ecosystem(proportions, Z)), [0, 1, 2, Inf])
diversity = supercommunityDiversity(RawRho(Ecosystem(proportions, Z)), 2)
```
