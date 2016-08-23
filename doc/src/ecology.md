In the **Diversity.Ecology** package, we replicate old ecological
diversity measures and generalised versions of them that relate to our
general measures of alpha, beta and gamma diversity at subcommunity
and ecosystem measures. The generalisations of the richness, Shannon
and Simpson are the only standard measures we are aware of whose
subcommunity components sum directly to the corresponding ecosystem
measure (although note that Simpson's index decreases for increased
diversity, so small components are more diverse).

#### Usage

Accessing the functionality in the package is simple:

```julia_skip
using Diversity.Ecology

community = [10. 20. 20.]'
diversity = simpson(community)

ecosystem = [2. 2. 0.; 0. 2. 2.]'
Z = eye(3)

jaccard(ecosystem)
generalisedjaccard(ecosystem, [0, 1, 2])
generalisedjaccard(ecosystem, [0, 1, 2], Z)
```

```@contents
```

```@autodocs
Modules = [Diversity.Ecology]
Private = false
```

```@index
```
