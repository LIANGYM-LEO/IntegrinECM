\name{score_integrin_ecm}
\alias{score_integrin_ecm}
\alias{score_integrin_ecm_files}
\alias{combine_integrin_ecm_results}
\alias{print.integrin_ecm_result}
\title{Run modular Integrin-ECM scoring}
\usage{
score_integrin_ecm(
  expression,
  metadata = NULL,
  modules = list(
    integrin = "nnls",
    collagen = "nnls",
    laminin = "geomean",
    other_ecm = "geomean"
  ),
  catalogs = NULL,
  missing_component_policy = "error",
  nnls_tolerance = 1e-10,
  output_digits = 2L,
  rank_ties_method = "average",
  include_zero_scores_in_ranking = TRUE,
  group_summary_function = "median"
)

score_integrin_ecm_files(
  expression_file,
  sample_metadata_file = NULL,
  gene_id = "gene",
  sample_id = "sample_id",
  group = "analysis_group",
  duplicate_rule = "error",
  ...
)

combine_integrin_ecm_results(...)
}
\arguments{
  \item{expression}{Non-negative gene-by-sample expression matrix.}
  \item{metadata}{Optional sample metadata.}
  \item{modules}{Named method-object selection list.}
  \item{catalogs}{Optional replacement definition list.}
  \item{missing_component_policy}{Missing-gene policy.}
  \item{nnls_tolerance}{NNLS numerical tolerance.}
  \item{output_digits}{Percentage-score decimal places.}
  \item{rank_ties_method}{Within-sample tie method.}
  \item{include_zero_scores_in_ranking}{Whether zero scores receive ranks.}
  \item{group_summary_function}{Group summary function.}
  \item{expression_file}{Expression CSV path.}
  \item{sample_metadata_file}{Optional metadata CSV path.}
  \item{gene_id}{Expression gene-column name.}
  \item{sample_id}{Metadata sample-column name.}
  \item{group}{Metadata group-column name.}
  \item{duplicate_rule}{Duplicate-gene rule.}
  \item{...}{For \code{score_integrin_ecm_files}, additional arguments passed to \code{score_integrin_ecm}. For \code{combine_integrin_ecm_results}, two or more result objects or one list of result objects.}
}
\description{Runs selected method-object modules, records reproducibility information, or combines compatible structured results.}
\value{An \code{integrin_ecm_result}.}
\details{When results are combined, sample metadata and scoring parameters must match. If expression-file checksums are available, they must also agree. Duplicate method-object modules are rejected.}
