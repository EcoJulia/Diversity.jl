The **Diversity** package provides functionality for measuring alpha,
beta and gamma diversity of subcommunities and ecosystems. It uses
diversity measures extended from those described in the arXiv paper
[arXiv:1404.6520 (q-bio.QM)](http://arxiv.org/abs/1404.6520),
*How to partition diversity*. The alpha, beta and gamma diversities in
the paper are supplemented by beta diversity measures related to
redundancy, DR and representativeness, DR̄, of ecosystems, and their associated
subcommunities (Dρ and Dρ̄, respectively).

# Usage

Accessing the functionality in the package is simple:

```julia_skip
using Diversity

# Load up ecosystem

diversities = Dᾱ(ecosystem, [0, 1, 2, Inf], Z)
diversities = Dγ(ecosystem, [0, 1, 2, Inf], Z)
```
