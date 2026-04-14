# ReMIXTURE VCF/BCF Distance-Matrix Change Report

Generated on: 2026-04-14

## 1. Scope of This Change

This update adds a new input path to `ReMIXTURE` so that a pairwise distance matrix can be built directly from genotype calls stored in a `VCF` or `BCF` file, instead of requiring the user to provide a precomputed `distance_matrix`.

The implementation was added primarily in `R/class.R` and includes:

- validation of `sample_table` and streaming parameters
- direct `GT` parsing and dosage conversion
- chunked streaming for large variant files
- pairwise distance accumulation across chunks
- optional local cache reuse for repeated test runs
- a plotting robustness fix for zero-overlap runs

As of the current working tree, the code changes affect 3 tracked source files, and this report is an additional documentation file.

## 2. Files Modified, Line Numbers, and What Changed

### 2.1 `DESCRIPTION`

File: `ReMIXTURE/DESCRIPTION`

Relevant lines:

- `30-34`

Change summary:

- Added `Rsamtools` and `VariantAnnotation` to `Suggests`.
- Adjusted the `rmarkdown` line so the new entries can be appended cleanly.

Current content around the modified block:

- `30`: `Suggests:`
- `31`: `knitr,`
- `32`: `rmarkdown,`
- `33`: `Rsamtools,`
- `34`: `VariantAnnotation`

Why it matters:

- `VariantAnnotation` is used for VCF streaming and `GT` extraction.
- `Rsamtools` is used for indexed BCF access.

### 2.2 `R/class.R`

File: `ReMIXTURE/R/class.R`

This is the main implementation file for the VCF/BCF feature. The change can be divided into 4 blocks.

#### A. New private VCF/BCF distance-matrix pipeline

Relevant block:

- `207-998`

New helper functions and start lines:

- `207` `validate_variant_sample_table()`
- `251` `validate_variant_cache_params()`
- `269` `build_variant_cache_metadata()`
- `289` `cache_metadata_matches()`
- `306` `maybe_load_cached_distance_matrix()`
- `334` `maybe_save_cached_distance_matrix()`
- `371` `validate_variant_stream_params()`
- `393` `detect_variant_file_type()`
- `407` `extract_header_ids()`
- `433` `normalise_variant_gt_chunk()`
- `486` `gt_to_dosage_matrix()`
- `523` `initialise_distance_accumulator()`
- `546` `accumulate_dosage_chunk()`
- `583` `finalize_distance_accumulator()`
- `637` `extract_bcf_seq_table()`
- `689` `build_bcf_chunk_ranges()`
- `713` `stream_vcf_distance_matrix()`
- `792` `stream_bcf_distance_matrix()`
- `913` `distance_matrix_from_vcf_inputs()`

What was added:

- A validation layer for `sample_table` to ensure:
  - the sample column exists
  - the region column exists
  - sample names are not missing
  - sample names are unique
  - sample names requested by the user exist in the VCF/BCF header
- A cache layer controlled by `cache_file` and `rebuild_cache`, keyed by:
  - normalized input file path
  - file size
  - modification time
  - ordered sample/region mapping
  - `min_variant_call_rate`
  - `min_pairwise_overlap`
- File-type detection so the code dispatches to either:
  - `stream_vcf_distance_matrix()`
  - `stream_bcf_distance_matrix()`
- `GT` normalization and dosage conversion rules:
  - `0/0 -> 0`
  - `0/1` or `1/0 -> 1`
  - `1/1 -> 2`
  - all other `GT` patterns are treated as missing for this workflow
- Variant-level filtering using `min_variant_call_rate`
- Pairwise overlap validation using `min_pairwise_overlap`
- Chunked VCF reading via `VariantAnnotation::VcfFile(..., yieldSize = chunk_variants)`
- Chunked BCF scanning via:
  - `Rsamtools::scanBcfHeader()`
  - contig-length parsing
  - genomic windows generated from `chunk_bp`
  - repeated `scanBcf()` calls on indexed windows
- Final pairwise distance assembly by accumulating mean absolute dosage differences over all valid shared sites

What this achieves:

- `ReMIXTURE` can now build a distance matrix directly from raw genotype files.
- The implementation scales better to larger inputs because it does not require loading the entire VCF/BCF into memory at once.
- Repeated debugging or plotting runs can reuse a cached matrix instead of re-streaming the variant file every time.

#### B. Extended `initialize()`

Relevant block:

- `1040-1147`

Key lines:

- `1044-1053`: new roxygen parameter documentation
- `1055-1067`: new arguments in `initialize()`
- `1077-1081`: mutual-exclusion checks between `distance_matrix` and `vcf_file`
- `1083-1099`: on-the-fly matrix construction from `vcf_file`

New arguments:

- `vcf_file`
- `sample_table`
- `sample_col`
- `region_col`
- `min_variant_call_rate`
- `min_pairwise_overlap`
- `chunk_variants`
- `chunk_bp`
- `cache_file`
- `rebuild_cache`

Behavioral change:

- Before this update, object construction depended on a precomputed `distance_matrix`.
- After this update, the user may provide either:
  - `distance_matrix`, or
  - `vcf_file` plus `sample_table`

Result:

- The package now supports direct initialization from raw genotype input:
  `ReMIXTURE$new(vcf_file = ..., sample_table = ..., region_table = ...)`

#### C. New public helper `distance_matrix_from_vcf()`

Relevant block:

- `1149-1192`

Key lines:

- `1156-1165`: roxygen documentation for the new public method
- `1168-1192`: method implementation

What was added:

- A public method `distance_matrix_from_vcf()` exposing the same main controls as the private streaming pipeline.

Why it matters:

- Users can now generate and inspect a distance matrix from a VCF/BCF file without constructing a full `ReMIXTURE` object first.
- This also makes the new feature easier to test in isolation.

#### D. Zero-overlap plotting fix in `plot_maps()`

Relevant block:

- `1670-1678`

What changed:

- The alpha scaling logic now checks whether the overlap-derived denominator is finite and positive before rescaling transparency.
- If the overlap matrix has no valid positive source value, alpha is set to `0.0` instead of producing `NaN`.

Why it matters:

- Runs with no between-region overlap now produce invisible lines instead of a plotting error caused by invalid alpha values.

### 2.3 `R/utility.R`

File: `ReMIXTURE/R/utility.R`

Relevant block:

- `119-126`

What changed:

- Added an explicit finite-value check for `setAlpha` inside `parseColChain()`.
- The new guard raises:
  - `"\'setAlpha\' values must be finite."`

Why it matters:

- This closes the plotting robustness loop by preventing hidden failures in color conversion when alpha carries `NA`, `Inf`, or `-Inf`.

## 3. Added Files

### 3.1 Added documentation file

- `ReMIXTURE/VCF_BCF_distance_matrix_change_report.md`

Purpose:

- Human-readable summary of the VCF/BCF feature change, dependencies, local test assets, reproducible scripts, and PR wording.

### 3.2 Added source files for the feature

- No new tracked R source files were added.
- The feature was implemented by extending existing files:
  - `DESCRIPTION`
  - `R/class.R`
  - `R/utility.R`

## 4. New Dependencies Required by This Feature

### 4.1 Explicitly declared in `DESCRIPTION`

- `VariantAnnotation`
- `Rsamtools`

### 4.2 Runtime namespaces directly referenced by the BCF code path

- `GenomicRanges`
- `IRanges`

Important note:

- `GenomicRanges::GRanges()` and `IRanges::IRanges()` are called directly in `R/class.R` at lines `856` and `858`.
- These packages are not currently declared in `DESCRIPTION`.
- In practice they are commonly installed together with Bioconductor tooling, but they are still runtime requirements for the BCF path.

### 4.3 Recommended installation command for local testing

```r
install.packages(c("remotes", "devtools", "data.table"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation", "Rsamtools", "GenomicRanges", "IRanges"))
```

## 5. Test Files Used in This Local Development Cycle

This section is based on files currently present on `D:/` in the local machine and on the interactive testing steps recorded during this development session.

### 5.1 Official upstream files used as test inputs

#### A. Full chr22 VCF from 1000 Genomes Phase 3

- Local file:
  - `D:/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`
- Local index:
  - `D:/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi`
- File size observed locally:
  - `214,453,750` bytes
- Role in testing:
  - main real-world VCF input for local end-to-end testing of `VCF -> distance_matrix -> ReMIXTURE`
- Upstream source:
  - 1000 Genomes Project Phase 3 release directory dated `20130502`
  - official release announcement: <https://www.internationalgenome.org/category/phase-3/>
- Exact release-path statement:
  - the exact filename used locally matches the official Phase 3 release naming convention under `release/20130502/`
  - the precise file URL is therefore reasonably inferred as:
    `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`

#### B. Official sample panel file

- Local file:
  - `D:/integrated_call_samples_v3.20130502.ALL.panel`
- File size observed locally:
  - `55,156` bytes
- Role in testing:
  - used to map sample IDs to population labels such as `CHB`, `GBR`, and `YRI`
- Upstream source:
  - official 1000 Genomes Phase 3 panel file
  - official path exposed through the project FTP tree:
    `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_call_samples_v3.20130502.ALL.panel`

### 5.2 Local derived files created for testing

These are not external published files. They were created locally from the official 1000 Genomes inputs above for faster or more focused testing.

#### A. 60-sample subset definition

- Local file:
  - `D:/samples.txt`
- Role:
  - sample ID list used to derive a small chr22 test VCF
- Observed content pattern:
  - 60 sample IDs total
  - 20 `CHB`
  - 20 `GBR`
  - 20 `YRI`

#### B. Small derived chr22 VCF for smoke tests

- Local file:
  - `D:/1kg_chr22_test.vcf.gz`
- Local index:
  - `D:/1kg_chr22_test.vcf.gz.tbi`
- File size observed locally:
  - `561,462` bytes
- Role:
  - small fast-running VCF used for smoke tests during debugging
- Provenance:
  - locally derived subset of the official full chr22 VCF above
  - not a separately published upstream dataset

#### C. Local sample/region helper tables

- `D:/sample_table.rds`
  - 3-population sample mapping used in the 60-sample test
- `D:/region_table.rds`
  - region centroid table used for plotting the 3-population test
- `D:/sample_table_all26.rds`
  - broader sample mapping for all 26 1000 Genomes populations
- `D:/region_table_all26.rds`
  - region centroid table for the all-26-population plotting test

Important note:

- `sample_table*.rds` files are local convenience products derived from the official panel file.
- `region_table*.rds` files are local helper metadata for plotting and are not official 1000 Genomes release files.

### 5.3 BCF test-file status in this local session

- No `.bcf` file is currently present in `D:/`.
- Therefore, the code path for `BCF` support was implemented and dependency-wired, but the local test evidence preserved on disk for this session is VCF-based rather than BCF-based.

This should be stated honestly in any PR or report.

## 6. Source Basis and Literature Citations

### 6.1 Dataset/source references for the local VCF test

1. 1000 Genomes Project Phase 3 release announcement  
   Link: <https://www.internationalgenome.org/category/phase-3/>

2. Official 1000 Genomes sample panel file  
   Link: <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_call_samples_v3.20130502.ALL.panel>

3. Official 1000 Genomes FTP release tree used for the chr22 VCF  
   Release directory announced for `20130502`; the locally used chr22 VCF filename matches this release naming scheme.

### 6.2 Paper citation for the 1000 Genomes data used in testing

Auton A, Brooks LD, Durbin RM, Garrison EP, Kang HM, Korbel JO, Marchini JL, McCarthy S, McVean GA, Abecasis GR, et al. (2015). *A global reference for human genetic variation*. Nature, 526, 68-74.  
DOI: `10.1038/nature15393`

Interpretation:

- This is the paper-level citation for the 1000 Genomes Phase 3 reference data underlying the local VCF and panel files used for testing.

### 6.3 Existing ReMIXTURE package citation context

Tripodi P, Rabanus-Wallace MT, et al. (2021). *Global range expansion history of pepper (Capsicum spp.) revealed by over 10,000 genebank accessions*. PNAS.  
DOI: `10.1073/pnas.2104315118`

Interpretation:

- This citation is relevant to the package's existing bundled example data and original concept history.
- It is not the literature source of the 1000 Genomes VCF used to test the new VCF/BCF import pathway.

## 7. Reproducible Local Test Scripts

The repository does not currently contain a committed automated test file for the new feature. The following scripts are a consolidated, reproducible version of the local testing workflow used in RStudio.

### 7.1 Install and load the local package

```r
install.packages(c("remotes", "devtools", "data.table"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation", "Rsamtools", "GenomicRanges", "IRanges"))

devtools::load_all("C:/Users/42006/Desktop/remixture/ReMIXTURE")
```

### 7.2 Main 60-sample chr22 VCF test

```r
library(data.table)

vcf_file <- "D:/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
sample_table <- readRDS("D:/sample_table.rds")
region_table <- readRDS("D:/region_table.rds")

rm_vcf <- ReMIXTURE$new(
  vcf_file = vcf_file,
  sample_table = sample_table,
  region_table = region_table,
  min_variant_call_rate = 0.95,
  min_pairwise_overlap = 1000L,
  chunk_variants = 10000L,
  cache_file = "D:/rm_chr22_cache.rds",
  rebuild_cache = FALSE
)

dim(rm_vcf$distance_matrix)
table(colnames(rm_vcf$distance_matrix))
range(rm_vcf$distance_matrix)

rm_vcf$plot_MDS()

rm_vcf$run(
  iterations = 10,
  subsample_proportions = 0.8
)

rm_vcf$plot_h_optimisation()
rm_vcf$plot_results_grid()
rm_vcf$plot_maps(
  run = 1,
  focalRegion = "GBR",
  width_max = 8
)
```

### 7.3 Optional fast smoke test using the local subset VCF

```r
vcf_file <- "D:/1kg_chr22_test.vcf.gz"
sample_table <- readRDS("D:/sample_table.rds")
region_table <- readRDS("D:/region_table.rds")

rm_small <- ReMIXTURE$new(
  vcf_file = vcf_file,
  sample_table = sample_table,
  region_table = region_table,
  min_variant_call_rate = 0.90,
  min_pairwise_overlap = 100L,
  chunk_variants = 5000L,
  cache_file = "D:/rm_chr22_small_cache.rds",
  rebuild_cache = FALSE
)

dim(rm_small$distance_matrix)
range(rm_small$distance_matrix)
```

### 7.4 Optional broader all-26-population test

```r
vcf_file <- "D:/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
sample_table <- readRDS("D:/sample_table_all26.rds")
region_table <- readRDS("D:/region_table_all26.rds")

rm_all26 <- ReMIXTURE$new(
  vcf_file = vcf_file,
  sample_table = sample_table,
  region_table = region_table,
  min_variant_call_rate = 0.95,
  min_pairwise_overlap = 1000L,
  chunk_variants = 10000L,
  cache_file = "D:/rm_chr22_all26_cache.rds",
  rebuild_cache = FALSE
)

dim(rm_all26$distance_matrix)
table(colnames(rm_all26$distance_matrix))
```

### 7.5 Optional script to regenerate the 60-sample `sample_table`

```r
library(data.table)

panel <- fread("D:/integrated_call_samples_v3.20130502.ALL.panel")
keep_samples <- fread("D:/samples.txt", header = FALSE)$V1

sample_table <- panel[
  sample %in% keep_samples,
  .(sample, region = pop)
]

saveRDS(sample_table, "D:/sample_table.rds")
```

## 8. Test Notes and Practical Cautions

1. `vcf_file` and `distance_matrix` are mutually exclusive. One must be provided, but not both.
2. If `vcf_file` is used, `sample_table` is required.
3. The sample column in `sample_table` must match the sample IDs in the VCF/BCF header exactly.
4. Only `GT` is used for the new workflow. Files without a valid `GT` field will fail.
5. The current implementation only interprets diploid biallelic dosage-like `GT` states:
   `0/0`, `0/1`, `1/0`, `1/1`.
6. Other genotype encodings are treated as missing in the current implementation.
7. `min_variant_call_rate` that is too strict can filter out nearly all sites.
8. `min_pairwise_overlap` that is too strict can reject the run if some sample pairs retain too few shared variants.
9. For repeated debugging runs, use `cache_file`; if code or inputs changed, set `rebuild_cache = TRUE`.
10. For Windows paths in R, use forward slashes such as `D:/file.vcf.gz`.
11. The plot can still look visually uninformative if only one or a few population pairs dominate the overlap signal. That is a data/result issue, not necessarily a parser failure.
12. `width_max` that is set too large can make map links look like thick bars and hide geographic detail.
13. The BCF path requires an indexed BCF and the runtime availability of `Rsamtools`, `GenomicRanges`, and `IRanges`.
14. The current working tree has not regenerated the Rd documentation after adding the new public arguments and method. Running `devtools::document()` is recommended before final packaging or PR submission.

## 9. Suggested English PR / Branch Description

### 9.1 Suggested PR title

`Add chunked VCF/BCF distance-matrix construction with cache support`

### 9.2 Suggested PR summary

```text
This change adds direct VCF/BCF input support to ReMIXTURE so that a pairwise distance matrix can be built from genotype calls during object initialization or through a dedicated helper method.

Key updates:
- add a chunked VCF/BCF streaming pipeline for GT-based distance-matrix construction
- support direct initialization from vcf_file + sample_table
- add a public distance_matrix_from_vcf() helper
- add optional cache_file / rebuild_cache controls for repeated local tests
- harden plotting against zero-overlap alpha-scaling failures

Validation:
- locally tested with the 1000 Genomes Phase 3 chr22 VCF and official sample panel
- verified on a 60-sample CHB/GBR/YRI workflow and on broader local sample-table variants
- current local on-disk evidence covers the VCF workflow; separate indexed BCF validation is still recommended
```

### 9.3 Short branch description

```text
Implements chunked VCF/BCF-to-distance-matrix support in ReMIXTURE, adds reusable local cache support for repeated debugging, and fixes zero-overlap plotting failures.
```

## 10. One-Sentence Executive Summary

This update extends `ReMIXTURE` from a distance-matrix-only workflow to a raw-genotype workflow by adding chunked VCF/BCF parsing, dosage-based distance construction, optional cache reuse, and a plotting robustness fix, with local end-to-end validation performed on real 1000 Genomes Phase 3 chr22 VCF data.
