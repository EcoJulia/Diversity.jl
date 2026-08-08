# SPDX-License-Identifier: BSD-2-Clause

module DiversityPopGenExt

import Diversity
using Diversity.API

using PopGen
using DataFrames

"""
    GeneticVCF

Genetic similarity type built from a `PopGen.PopData` object (e.g. read from a
VCF file). Each sample is a type; similarity is derived from pairwise genetic
distances between samples.
"""
struct GeneticVCF{PopDataType} <: Diversity.AbstractGenetic
    dat::PopDataType
    names::Vector{String}
    ntypes::Int64
    Zmatrix::Matrix{Float64}
end

# The reference allele for a locus: the smallest allele code present. Manhattan
# distance on allele dosage is invariant to this choice, so any consistent pick
# reproduces rdiversity's gen2dist(biallelic = TRUE).
function _reference_allele(column)
    ref = nothing
    for genotype in column
        ismissing(genotype) && continue
        for allele in genotype
            ref = isnothing(ref) ? allele : min(ref, allele)
        end
    end
    return ref
end

# Dosage = number of non-reference alleles in a genotype. Missing genotypes are
# treated as homozygous reference (dosage 0), matching gen2dist's recoding of
# missing data as no mutation.
_dosage(::Missing, ref) = 0
_dosage(genotype, ref) = isnothing(ref) ? 0 : count(!=(ref), genotype)

# Biallelic Manhattan distance between samples from a (samples × loci) genotype
# matrix, matching rdiversity gen2dist(vcf, biallelic = TRUE).
function _distancematrix(::Val{:manhattan}, lm::AbstractMatrix)
    nsamples, nloci = size(lm)
    dist = zeros(Float64, nsamples, nsamples)
    for locus in 1:nloci
        column = view(lm, :, locus)
        ref = _reference_allele(column)
        doses = [_dosage(column[s], ref) for s in 1:nsamples]
        for a in 1:nsamples, b in (a + 1):nsamples
            d = abs(doses[a] - doses[b])
            dist[a, b] += d
            dist[b, a] += d
        end
    end
    return dist
end

# Per-locus genotype Hamming distance between samples: the number of loci at
# which two samples carry a different (unordered) genotype.
function _distancematrix(::Val{:hamming}, lm::AbstractMatrix)
    nsamples, nloci = size(lm)
    dist = zeros(Float64, nsamples, nsamples)
    for a in 1:nsamples, b in (a + 1):nsamples
        d = 0
        for locus in 1:nloci
            ga, gb = lm[a, locus], lm[b, locus]
            differ = if ismissing(ga) || ismissing(gb)
                !(ismissing(ga) && ismissing(gb))
            else
                sort!(collect(ga)) != sort!(collect(gb))
            end
            d += differ
        end
        dist[a, b] = dist[b, a] = d
    end
    return dist
end

function _distancematrix(::Val{D}, ::AbstractMatrix) where {D}
    return throw(ArgumentError("unknown genetic distance :$D (try :manhattan or :hamming)"))
end

"""
    GeneticType(dat::PopData; distance = :manhattan, names = samplenames(dat),
                transform = :linear, k = 1, normalise = true)

Construct a `GeneticVCF` similarity type from a `PopGen.PopData` object. Samples
are the types. `distance` selects the pairwise distance (`:manhattan`, matching
rdiversity's `gen2dist(vcf, biallelic = TRUE)`, or `:hamming`), and `transform`,
`k` and `normalise` control the distance-to-similarity conversion (see
rdiversity's `dist2sim`).
"""
function Diversity.GeneticType(dat::PopData;
                               distance::Symbol = :manhattan,
                               names::AbstractVector = samplenames(dat),
                               transform::Symbol = :linear,
                               k::Real = 1, normalise::Bool = true)
    lm = locimatrix(dat)
    dist = _distancematrix(Val(distance), lm)
    Zmatrix = Diversity._dist2sim(dist; transform = transform, k = k,
                                  normalise = normalise, max_d = maximum(dist))
    return GeneticVCF{typeof(dat)}(dat, String.(names), size(lm, 1), Zmatrix)
end

# Reconstruct a VCF genotype string from a PopGen genotype tuple. PopGen stores
# alleles as 1-based codes (ref = 1, alt = 2, …), so mapping code k -> k-1 gives
# the original 0-based VCF genotype (e.g. (1, 2) -> "0|1"). Missing genotypes
# become ".|.", which gen2dist recodes as no mutation.
_genotype_string(::Missing) = ".|."
_genotype_string(genotype) = join((allele - 1 for allele in genotype), "|")

"""
    vcf_dataframe(dat::PopData)

Convert a `PopGen.PopData` object into a `DataFrame` laid out like the body of a
VCF file: a `FORMAT` column followed by one genotype column per sample (rows are
loci). This is the structure consumed by rdiversity's `gen2dist()`, so the same
`PopData` read into Julia can be handed to R (e.g. `@rput vcf_dataframe(pd)`) to
drive both Julia and R genetic diversity calculations from a single source.
"""
function Diversity.vcf_dataframe(dat::PopData)
    lm = locimatrix(dat)                       # samples × loci genotype tuples
    samples = String.(samplenames(dat))
    nloci = size(lm, 2)
    df = DataFrame(CHROM = "1", POS = 1:nloci, ID = ".", REF = "A", ALT = "T",
                   QUAL = ".", FILTER = ".", INFO = ".", FORMAT = "GT")
    for (s, name) in enumerate(samples)
        df[!, name] = [_genotype_string(lm[s, locus]) for locus in 1:nloci]
    end
    return df
end

import Diversity.API: _getdiversityname
_getdiversityname(::GeneticVCF) = "Genetic (VCF)"

end
