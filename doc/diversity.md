The **Diversity** package provides basic numbers-equivalent diversity
measures (described in
[Hill, 1973](http://www.jstor.org/stable/1934352)),
similarity-sensitive diversity measures (generalised from Hill, and
described in
[Leinster and Cobbold, 2012](http://www.esajournals.org/doi/abs/10.1890/10-2402.1)),
and related alpha, beta and gamma diversity measures at the level of
the supercommunity and its component subcommunities (generalised in
turn from Leinster and Cobbold, and described in
[Reeve et al, 2014](http://arxiv.org/abs/1404.6520)). The functions
exist both with unicode names preceded by D (e.g. Dᾱ()), and with
matching ascii names (e.g. subcommunityalphabar()). We also provide a
general function for extract any diversity measure for a series of
subcommunity relative abundances. The full documentation can be found
[here](http://diversityjl.readthedocs.org/en/stable/diversity/).

#### Usage

Accessing the functionality in the package is simple:

```julia_skip
using Diversity

# Load up ecosystem

diversities = Dᾱ(ecosystem, [0, 1, 2, Inf], Z)
diversities = Dγ(ecosystem, [0, 1, 2, Inf], Z)
```
