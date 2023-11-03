# Diversity.API

The **Diversity.API** submodule provides the API that must be extended
for new `AbstractTypes`, `AbstractPartition` and
`AbstractMetacommunity` subtypes.

## Usage

Providing additional code to extend the functionality of the system is simple:

```julia
using Diversity
using Phylo
import Diversity.API: _counttypes

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
Modules = [Diversity.API, Diversity.powermean, Diversity._getmeta]
Private = false
```

```@index
```
