\name{read_expression_matrix}
\alias{read_expression_matrix}
\alias{read_sample_metadata}
\alias{write_score_workbook}
\title{Input and workbook-output helpers}
\usage{
read_expression_matrix(path, gene_id = "gene", duplicate_rule = "error")

read_sample_metadata(path, sample_id = "sample_id", group = "analysis_group")

write_score_workbook(
  result,
  path,
  optional_sheets = character(0),
  font_name = NULL
)
}
\arguments{
  \item{path}{Input CSV path or output XLSX path, depending on the function.}
  \item{gene_id}{Expression gene-column name.}
  \item{duplicate_rule}{Duplicate-gene aggregation rule.}
  \item{sample_id}{Metadata sample-column name.}
  \item{group}{Metadata group-column name.}
  \item{result}{An \code{integrin_ecm_result}.}
  \item{optional_sheets}{Optional technical workbook sheets.}
  \item{font_name}{Deprecated compatibility argument. Font post-processing is no longer applied because rebuilding XLSX files was not reliable across operating systems.}
}
\description{Reads canonical CSV inputs and writes compact or expanded result workbooks.}
\value{An expression matrix, metadata data frame, or invisible sheet-name vector.}
