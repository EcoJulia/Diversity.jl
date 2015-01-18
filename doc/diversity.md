The **Diversity** package provides functionality for measuring alpha,
beta and gamma diversity of subcommunities and ecosystems. It uses
diversity measures extended from those described in the arXiv paper
[arXiv:1404.6520 (q-bio.QM)](http://arxiv.org/abs/1404.6520),
*How to partition diversity*. The alpha, beta and gamma diversities in
the paper are supplemented by beta diversity measures related to
evenness, DE, and redundancy, DR, of ecosystems, and their associated
subcommunities (Dϵ and Dρ, respectively).

# Usage

Accessing the functionality in the package is simple:

```julia
using Diversity
...
diversities = Dᾱ(proportions, [0, 1, 2, Inf], Z)
diversities = Dγ(proportions, [0, 1, 2, Inf], Z)
```
