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
  "IntegrinECM_0.1.0.tar.gz",
  repos = NULL,
  type = "source"
)
library(IntegrinECM)
```

## Default four-module analysis

```r
library(IntegrinECM)

E <- read_expression_matrix("expression.csv", gene_id = "gene")
metadata <- read_sample_metadata(
  "sample_metadata.csv",
  sample_id = "sample_id",
  group = "analysis_group"
)

result <- score_integrin_ecm(E, metadata)
print(result)

write_score_workbook(
  result,
  "Integrin_ECM_Percent_Scores.xlsx"
)
```

Default method-object combinations are:

- `integrin_nnls`
- `collagen_nnls`
- `laminin_geomean`
- `other_ecm_geomean`

Each percentage-score row sums to exactly `100.00` after balanced rounding.

## Input data

The package does not contain patient-level expression or clinical data. Users
provide a gene-by-sample expression CSV and, optionally, a sample metadata CSV.
If metadata are omitted, all samples are assigned to `All samples`.

Bundled reference resources define the supported Integrin, Collagen, Laminin,
and Other-ECM structures and the curated 284-record Integrin--ECM interaction
table. Collagen score-column labels display the chain stoichiometry, for example
`Collagen I [(COL1A1)2(COL1A2)1]`.

## Selective execution

```r
collagen_only <- score_integrin_ecm(
  E,
  metadata,
  modules = list(collagen = "nnls")
)

collagen_both <- score_integrin_ecm(
  E,
  metadata,
  modules = list(collagen = c("nnls", "geomean"))
)
```

## Custom structures

`score_structures()` accepts any compatible definition containing:

- `catalog`: `structure_id`, `structure_label`;
- `components`: `structure_id`, `structure_label`, `gene`, `coefficient`.

This keeps both methods general rather than restricted to Integrin or ECM.
