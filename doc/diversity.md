The **Diversity** package provides functionality for measuring alpha,
beta and gamma diversity of subcommunities and ecosystems. It uses
diversity measures extended from those described in the arXiv paper
[arXiv:1404.6520 (q-bio.QM)](http://arxiv.org/abs/1404.6520),
*How to partition diversity*.

# Usage

Accessing the functionality in the package is simple:

```julia
using Diversity
...
diversities = Diversity.ᾱ(proportions, [0, 1, 2, Inf], Z)
diversities = gamma(proportions, [0, 1, 2, Inf], Z)
```
