# Replace these paths with the user's own files.
expression_file <- "expression.csv"
sample_metadata_file <- "sample_metadata.csv"
output_file <- "Integrin_ECM_Percent_Scores.xlsx"

library(IntegrinECM)

E <- read_expression_matrix(expression_file, gene_id = "gene")
metadata <- read_sample_metadata(
  sample_metadata_file,
  sample_id = "sample_id",
  group = "analysis_group"
)

result <- score_integrin_ecm(
  expression = E,
  metadata = metadata,
  modules = list(
    integrin = "nnls",
    collagen = "nnls",
    laminin = "geomean",
    other_ecm = "geomean"
  )
)

write_score_workbook(result, output_file)
