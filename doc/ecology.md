In the **Diversity.Ecology** package, we replicate old ecological
diversity measures and generalised versions of them that relate to our
general measures of alpha, beta and gamma diversity at subcommunity
and ecosystem measures. The generalisations of the richness, Shannon
and Simpson are the only standard measures we are aware of whose
subcommunity components sum directly to the corresponding ecosystem
measure (although note that Simpson's index decreases for increased
diversity, so small components are more diverse).

# Usage

Accessing the functionality in the package is simple:

```julia
using Diversity.Ecology
...
diversities = simpson(proportions)
diversities = generalisedjaccard(proportions, [0, 1, 2], Z)
```
