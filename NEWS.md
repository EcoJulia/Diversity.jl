# NEWS

- v0.6.2
  - Make measurement dramatically faster and less memory-hungry fixing getASCIIName()
  - Fix typematch() rejecting float types that are not direct subtypes of AbstractFloat
  - Build the result DataFrame from whole columns rather than one single-row DataFrame per element
  - Allow the caller to name a Tables.jl sink as an optional first argument, as CSV.read() does, so
    results can go to CSV or Arrow without a DataFrame in between; DataFrame remains the default
  - Hold each measure's individual diversities as a rule for computing an element rather than as an
    array to remove the largest allocation in a measurement, at the cost of about a sixth of the
    analysis time, since each element is recomputed rather than re-read every time it is used
  - Cache the subcommunity weights and the metacommunity ordinariness on a Metacommunity, as the
    ordinariness itself already was
  - Complex diversity calculations can create huge outputs, so there are now generators for the
    repeated columns since enormous storage costs until they are materialised (which isn't needed
    if writing to disk)
  - Recognise an empty subcommunity from its cached weight instead of by scanning its abundances,
    so a metacommunity that is mostly empty costs very nearly what it would if those subcommunities
    were not there at all
  - Add a manual page on large metacommunities, covering what actually scales with the data, using
    diversity() to ask for everything at once, sending results straight to a sink rather than
    through a DataFrame, and why empty subcommunities are nearly free
  - Document the notation convention on the framework page: a symbol's alphabet and case give its
    level, lowercase Roman for individual values, lowercase Greek for subcommunity measures and
    capital Roman for metacommunity ones
- v0.6.1
  - Implement EcoBase's view() for metacommunities, returning a lazy SubAssemblage that aliases the
    parent's abundances; this also makes EcoBase's cooccurring() and SpatialEcology's groupsites()
    and groupspecies() work on a metacommunity, which were previously a MethodError
  - Add _subsettypes() and _subsetpartition() to the extension contract, as optional functions with
    defaults, so a new type can say how it is subset
  - Fix showing a metacommunity using the singular where it needed the plural, and name the units
    through EcoBase's thingkind()/placekind() hooks, so output reads "3 species in 4 subcommunities"
    rather than "3 thing in 4 place", and says "branches" for a phylogeny
- v0.6.0
  - Document the diversity framework and its origins, with a new page in the manual
  - Add a page on building metacommunities, and one for users coming from R's vegan
  - Document the options to gower(), and how they relate to vegan's gower and altGower
  - Document the fourteen norm_/raw_ diversity wrapper functions, which had no docstrings
  - Document what each optional function of the extension contract does if you do not implement it
  - Document the exported AbstractPhyloTypes and PhyloBranches, which had no docstrings
  - Build the manual with Documenter in the test suite, so a dangling cross-reference fails there
    rather than only in the documentation workflow
  - Rename the getFullName() descriptions of the measures to match those in the paper
  - Fix the NormalisedRho docstring, which named it redundancy rather than representativeness
  - Fix docstrings giving the maximum of raw rho and normalised beta as the number of
    subcommunities, rather than the effective number of them
  - Add Faith's PD to Diversity.Ecology as faith_pd() and generalisedfaith_pd()
  - Add in genetic diversity, matching the gen2dist() and dist2sim() pipeline in boydorr/rdiversity
  - Provide it through two new extensions, BioSequences for aligned sequences and PopGen for VCF data
  - Cross-validate genetic diversity against rdiversity for every VCF file in test/data
  - Fix showing a metacommunity whose types are not named with strings, such as Metacommunity(pop, Z)
  - Fix the plot recipes, which errored on any input, and add tests for them
  - Fix Metacommunity(assemblage) for an assemblage whose types carry similarity, if any exist!
  - Fix a stack overflow when a type subtyping AbstractMetacommunity, AbstractPartition or
    AbstractTypes implements neither the Diversity API nor the EcoBase interface, which now
    reports the missing method instead
  - Report an AbstractTypes that says it has similarity but never implements _calcsimilarity(),
    rather than silently measuring it with an identity matrix; types that declare
    _hassimilarity() false still get that matrix as their default
  - Allow Metacommunity(counts, Z) with integer counts, as every other constructor does
  - Move StringDistances into core dependencies
  - Remove pre-v1.9 Requires syntax
  - Update compat for SpatialEcology and JuliaFormatter
- v0.5.16
  - Fix problem with R vegan validation
  - Update hygiene tests
  - Add DOI for arXiv paper
- v0.5.15
  - Use ResearchSoftwareMetadata package
  - Extensions of crosswalk, bugfixes
  - Add in testing that metadata is up-to-date
- v0.5.14
  - Add in crosswalk between Project.toml and codemeta.json, .zenodo.json, LICENSE and Julia file headers
- v0.5.13
  - Fix CI for arm64
  - Add in codemeta.json to allow metadata to follow RSMD standards
- v0.5.12
  - Move Phylo structs into extension
- v0.5.11
  - Update Phylo compat
- v0.5.10
  - Use extensions for Julia 1.9+
  - Allow AxisArrays to name types through extension
  - Improve testing
- v0.5.9
  - Introduce Gower
- v0.5.8
  - Add in Pielou diversity measures
  - Add docs for PRs
  - Generalise Jaccard with similarity
- v0.5.7
  - Update compat and some doc fixes
- v0.5.6
  - Fix phylogenetics
- v0.5.5
  - Update phylogenetics code and docs
- v0.5.4
  - Move to EcoJulia
  - Update compats
- v0.5.3
  - Update docs and automation
- v0.5.2
  - Update docs
  - Update compats
- v0.5.1
  - Improve dependencies
- v0.5.0
  - Drop Julia pre-v1.0 compatibility
  - Move to Phylo v0.4.0 API
- v0.4.5
  - Full EcoJulia default implementation to allow cross-compatibility with other AbstractAssemblage-derived types
- v0.4.4
  - Major reworking, especially of iterators to speed up code and run with Julia v1.0
  - Some intermediate work to allow interface to EcoJulia and new PhyloSets
- v0.4.3
  - Minor fixes to validate with updates to Julia and rdiversity R package
  - Fix some more Julia 0.7 problems
- v0.4.2
  - Work around problems with PackageEvaluator.jl
  - Julia 0.7 iterator fixes
- v0.4.1
  - Fix Nullable{T} to Union{T, Missing}
  - Other minor updates for Julia 0.7
- v0.4.0
  - Add in phylogenetic diversity
  - Create a formal API (in `API.jl`) for extending to new types of diversity
  - Extract interface into `Interface.jl`
  - Julia v0.6 and nightly compliant, drop support for Julia v0.5
  - Removes deprecated syntax from Diversity v0.2.x
- v0.3.1
  - Include validation against boydorr/rdiversity R package
- v0.3.0
  - Update input interface and deprecate old format
  - Update output format to use DataFrames
  - Remove deprecations for Julia v0.6 and drop support for Julia v0.4
