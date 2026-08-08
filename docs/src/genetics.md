# Genetic diversity

Genetic diversity is provided by two lightweight extensions, so you only load
what your data needs:

- **sequences** — load `BioSequences` (`using Diversity, BioSequences`) to build a
  `GeneticType` from a vector of aligned `BioSequence`s;
- **VCF** — load `PopGen` (`using Diversity, PopGen`) to build a `GeneticType`
  from a `PopGen.PopData` object read from a VCF file.

Similarity between types is derived from pairwise genetic distances, mirroring
the `gen2dist()` / `dist2sim()` pipeline in the R package
[rdiversity](https://github.com/boydorr/rdiversity).

## Usage

Using the functionality in the package is simple:

- Create genetic data, either a vector of aligned `BioSequence`s or a
  `PopGen.PopData` object (e.g. read from a VCF file with `PopGen.vcf`)
- Create a `GeneticType` (an `AbstractTypes` subtype) from it
- Create a `Metacommunity` from that
- Calculate diversity!

### Sequences

```@repl sequence
using Diversity, BioSequences
seqs = [dna"ACGTACGT", dna"ACGAACGT", dna"TTTTTTTT"]
gt = GeneticType(seqs; names = ["a", "b", "c"])
calcsimilarity(gt, 1.0)
metagen = Metacommunity([0.3, 0.3, 0.4], gt)
meta_gamma(metagen, 0)
```

Sequences `a` and `b` differ at one site out of eight and are correspondingly
similar; both are far from `c`, which shares no site with either. The
metacommunity therefore holds rather less than three types' worth of diversity.

### VCF

The example below reads the small biallelic VCF that ships with the package, and
builds a similarity matrix using biallelic Manhattan distances — matching
rdiversity's `gen2dist(vcf, biallelic = TRUE)`:

```@repl vcf
using Diversity, PopGen
pd = PopGen.vcf(Diversity.path("data", "biallelic.vcf"); silent = true)
gt = GeneticType(pd; distance = :manhattan)
gettypenames(gt, true)
calcsimilarity(gt, 1.0)
metagen = Metacommunity([0.2 0.1; 0.1 0.3; 0.2 0.1], gt)
meta_gamma(metagen, 0)
```

Here the *samples* are the types, so a metacommunity is a matrix of sample
abundances with one column per subcommunity.

The `distance` (`:manhattan` or `:hamming`), `transform` (`:linear` or
`:exponential`), `k` and `normalise` keyword arguments control how the pairwise
distances and the resulting similarity matrix are calculated.

`vcf_dataframe(pd)` converts a `PopData` back into the VCF-body layout that
rdiversity's `gen2dist()` consumes, so the same data can drive both the Julia
and the R calculation — which is how the two are cross-validated against each
other in `test/run_rcall.jl`.

```@contents
```

```@index
```
