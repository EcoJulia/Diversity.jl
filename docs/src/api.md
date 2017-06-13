The **Diversity.API** submodule provides the API that must be extended
for new `AbstractTypes`, `AbstractPartition` and
`AbstractMetacommunity` subtypes.

#### Usage

Providing additional code to extend the functionality of the system is simple:

```
using Diversity.Phylogenetics
importall Diversity.API

function _counttypes(phy::Phylogeny)
    return phy.nancestral
end
```

extends `Diversity.API._counttypes()` (and therefore the directly
accessible `counttypes()` interface) to handle the `Phylogeny` subtype
of `AbstractTypes`.

```@contents
```

```@autodocs
Modules = [Diversity.API]
Private = false
```

```@index
```
