# Diversity.Hill

[Hill numbers](http://www.jstor.org/stable/1934352) are found in the
**Diversity.Hill** submodule.

## Usage

Accessing the main functionality in the package is simple:

```@repl hill
using Diversity.Hill
community = [10, 20, 20, 0, 3]
community = community ./ sum(community)
hillnumber(community, [0, 1, 2])
```

A Hill number of order `q` is the diversity of a community in which every type
is completely distinct from every other. Order `0` counts the types present -
here `4`, since one of the five has no individuals - and higher orders weight
the commoner types more heavily, so the diversity falls as `q` rises.

```@contents
```

```@autodocs
Modules = [Diversity.Hill]
Private = false
```

```@index
```
