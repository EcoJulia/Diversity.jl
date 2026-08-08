# Diversity.Jost

Lou Jost's
[diversity](http://dx.doi.org/10.1111/j.2006.0030-1299.14714.x)
[measures](http://www.esajournals.org/doi/abs/10.1890/06-1736.1) are
found in the **Diversity.Jost** submodule.

## Usage

Accessing the main functionality in the package is simple:

```@repl jost
using Diversity.Jost
ecosystem = [2 2 0; 0 2 2]'
ecosystem = ecosystem ./ sum(ecosystem)
jostbeta(ecosystem, [0, 1, 2])
jostalpha(ecosystem, [0, 1, 2])
```

Jost's beta is the naive gamma diversity divided by Jost's alpha, and his alpha
is in turn the raw alpha diversity divided by the naive-community beta. We
believe our own [`NormalisedBeta`](@ref) has better properties — in particular it
is invariant under shattering, which Jost's alpha and beta are not, as
[the framework](framework.md) explains — but these are provided for comparison.

```@contents
```

```@autodocs
Modules = [Diversity.Jost]
Private = false
```

```@index
```
