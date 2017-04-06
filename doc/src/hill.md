[Hill numbers](http://www.jstor.org/stable/1934352) are found in the
**Diversity.Hill** package.

#### Usage

Accessing the main functionality in the package is simple:

```julia
using Diversity.Hill

# Load community to study

community = [10, 20, 20, 0, 3];
community /= sum(community);
diversities = hillnumber(community, [0, 1, 2])
```

```@contents
```

```@autodocs
Modules = [Diversity.Hill]
Private = false
```

```@index
```
