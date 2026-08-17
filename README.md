# IntegrinECM

`IntegrinECM` is an R implementation of a structure-aware scoring algorithm for
Integrin and extracellular-matrix assemblies. It accepts non-negative,
linear-scale gene-expression values and predefined component stoichiometry.

## Mathematical functions

The public functions retain the manuscript notation:

- `score_nnls_sample(A, E_s)` estimates `x_hat_s` and residual `epsilon_s`.
- `score_geomean_structure(t_s_i, v_i_k)` returns `u_s_k` and `q_s_k`.
- `normalize_structure_scores(h_s_k)` returns `p_s_k`, `z_s_k`, and ranks.

No pseudocount is used. If any required component has zero expression, the
geometric-mean score for that structure is zero.

## Installation from GitHub

```r
install.packages("remotes")
remotes::install_github("LIANGYM-LEO/IntegrinECM")
library(IntegrinECM)
```

## Installation from a local source archive

Download the release archive and install it using its full path:

```r
install.packages(
  "IntegrinECM_0.1.2.tar.gz",
  repos = NULL,
  type = "source"
)
library(IntegrinECM)
```

## Default four-module analysis

```r
library(IntegrinECM)

result <- score_integrin_ecm_files(
  expression_file = "expression.csv",
  sample_metadata_file = "sample_metadata.csv",
  gene_id = "gene",
  sample_id = "sample_id",
  group = "analysis_group"
)
print(result)

write_score_workbook(
  result,
  "Integrin_ECM_Percent_Scores.xlsx"
)
```

Workbook output is written directly with `writexl` for portability across
macOS, Windows, and Linux. The package does not unpack and rebuild XLSX files
to force a system-specific font.

Default method-object combinations are:

- `integrin_nnls`
- `collagen_nnls`
- `laminin_geomean`
- `other_ecm_geomean`

Within each method-object module, the displayed percentages sum to exactly
`100.00` for every sample after balanced rounding.

## Input data

The package does not contain patient-level expression or clinical data. Users
provide a gene-by-sample expression CSV and, optionally, a sample metadata CSV.
If metadata are omitted, all samples are assigned to `All samples`.

Bundled reference resources define the supported Integrin, Collagen, Laminin,
and Other-ECM structures and the curated 284-record Integrin--ECM interaction
table. Collagen score-column labels display the chain stoichiometry, for example
`Collagen I [(COL1A1)2(COL1A2)1]`.

## Bundled 10-sample example

Two small CSV files provide a directly runnable input example. They contain a
reproducible random subset of 10 TCGA-BRCA solid-tissue Normal samples
derived from gene-expression data obtained through the NCI Genomic Data
Commons. The expression columns and metadata sample identifiers are stored in
the same order. These data are included only to demonstrate the software
workflow and are not intended as an analytical cohort.

```r
expression_file <- system.file(
  "extdata", "example_expression_10_samples.csv",
  package = "IntegrinECM"
)
metadata_file <- system.file(
  "extdata", "example_sample_metadata_10_samples.csv",
  package = "IntegrinECM"
)

example_result <- score_integrin_ecm_files(
  expression_file,
  metadata_file,
  gene_id = "gene",
  sample_id = "sample_id",
  group = "analysis_group"
)

write_score_workbook(
  example_result,
  "IntegrinECM_10_sample_example.xlsx"
)
```

The same workflow is available as
`system.file("examples", "run_bundled_example.R", package = "IntegrinECM")`.

## Reference resources and reproducibility

The scoring workflow uses the bundled 284-record Integrin-ECM interaction
resource together with fixed definitions for Integrin, Collagen, Laminin, and
Other ECM structures. These resources are available from R whenever you want
to inspect the candidate structures or their component stoichiometry:

```r
interaction_reference <- load_interaction_reference()
integrin_reference <- load_structure_catalog("integrin")
collagen_reference <- load_structure_catalog("collagen")
laminin_reference <- load_structure_catalog("laminin")
other_ecm_reference <- load_structure_catalog("other_ecm")

nrow(interaction_reference)       # 284
head(collagen_reference$catalog)
head(collagen_reference$components)
```

To review or share the complete reference collection, write it to a separate
workbook; no expression data or scoring run is required:

```r
write_reference_workbook("IntegrinECM_reference_resources.xlsx")
```

For a file-based analysis, `score_integrin_ecm_files()` keeps the normalized
input paths and MD5 checksums with the result. The same object also stores the
effective parameters, selected modules, gene coverage, completion time, and
run status, making it easier to trace how a workbook was produced:

```r
result$input_audit
result$parameters
result$modules
result$run_status
result$completed_at
result$results$collagen_nnls$coverage
```

If you prefer to keep the full result in R rather than expand it into many
worksheets, save it with `saveRDS(result, "IntegrinECM_result.rds")` and restore
it later with `readRDS()`.

## Workbook output levels

The default workbook is deliberately compact and is intended for routine
inspection. More detailed calculations and reference tables can be added only
when they are needed:

```r
# Add raw scores, fitted values, diagnostics, coverage, and other intermediates
write_score_workbook(
  result,
  "IntegrinECM_intermediates.xlsx",
  optional_sheets = "intermediates"
)

# Add the references used in this run together with the 284 interaction records
write_score_workbook(
  result,
  "IntegrinECM_references.xlsx",
  optional_sheets = "references"
)

# Export the complete technical record
write_score_workbook(
  result,
  "IntegrinECM_complete.xlsx",
  optional_sheets = "all"
)
```

For finer control, individual sheets can be requested with `raw_long`,
`summary`, `diagnostics`, `fitted`, `coverage`, `composition`,
`geomean_intermediates`, `audit`, `annotations`, `components`, or
`interaction_reference`.

## Selective execution

```r
expression <- read_expression_matrix("expression.csv", gene_id = "gene")
metadata <- read_sample_metadata("sample_metadata.csv")

integrin_only <- score_integrin_ecm(
  expression,
  metadata,
  modules = list(integrin = "nnls")
)

collagen_only <- score_integrin_ecm(
  expression,
  metadata,
  modules = list(collagen = "nnls")
)

collagen_both <- score_integrin_ecm(
  expression,
  metadata,
  modules = list(collagen = c("nnls", "geomean"))
)
```

Integrin and Collagen can each use NNLS, the geometric mean, or both. Laminin
and Other ECM use the geometric mean in the current implementation:

```r
sensitivity_result <- score_integrin_ecm(
  expression,
  metadata,
  modules = list(
    integrin = c("nnls", "geomean"),
    collagen = c("nnls", "geomean"),
    laminin = "geomean",
    other_ecm = "geomean"
  )
)
```

Results produced in separate calls can also be brought together. Before doing
so, the package checks that the sample metadata and all scoring and
normalization settings agree; file-based runs must also have matching
expression checksums. It rejects duplicate method-object modules instead of
silently replacing one result with another:

```r
combined_result <- combine_integrin_ecm_results(
  integrin_run = integrin_only,
  collagen_run = collagen_only
)
```

## Custom structures

`score_structures()` accepts any compatible definition containing:

- `catalog`: `structure_id`, `structure_label`;
- `components`: `structure_id`, `structure_label`, `gene`, `coefficient`.

This keeps both methods general rather than restricted to Integrin or ECM.
