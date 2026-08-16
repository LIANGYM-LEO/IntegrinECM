\name{score_structures}
\alias{score_structures}
\alias{summarize_structure_scores}
\alias{print.integrin_ecm_module_result}
\title{Score one method-object module}
\usage{
score_structures(
  expression,
  definition,
  method = c("nnls", "geomean"),
  metadata = NULL,
  object = NULL,
  missing_component_policy = "error",
  nnls_tolerance = 1e-10,
  output_digits = 2L,
  rank_ties_method = "average",
  include_zero_scores_in_ranking = TRUE,
  group_summary_function = "median"
)

summarize_structure_scores(
  scores,
  summary_function = "median",
  ties_method = "average"
)
}
\arguments{
  \item{expression}{Non-negative gene-by-sample expression matrix.}
  \item{definition}{Compatible structure definition.}
  \item{method}{\code{"nnls"} or \code{"geomean"}.}
  \item{metadata}{Optional sample metadata.}
  \item{object}{Object label.}
  \item{missing_component_policy}{\code{"error"} or \code{"omit_structure"}.}
  \item{nnls_tolerance}{NNLS numerical tolerance.}
  \item{output_digits}{Percentage-score decimal places.}
  \item{rank_ties_method}{Within-sample tie method.}
  \item{include_zero_scores_in_ranking}{Whether zero scores receive ranks.}
  \item{group_summary_function}{\code{"median"} or \code{"mean"}.}
  \item{scores}{Long-format sample-level score table.}
  \item{summary_function}{Group summary function.}
  \item{ties_method}{Group-rank tie method.}
}
\description{
Runs one scoring method for one structure catalog, then normalizes and ranks
each sample independently.
}
\value{A module result or group-summary data frame.}
