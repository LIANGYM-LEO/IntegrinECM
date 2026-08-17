# Run the bundled, strictly matched 10-sample example.
library(IntegrinECM)

expression_file <- system.file(
  "extdata", "example_expression_10_samples.csv",
  package = "IntegrinECM"
)
sample_metadata_file <- system.file(
  "extdata", "example_sample_metadata_10_samples.csv",
  package = "IntegrinECM"
)
output_file <- file.path(getwd(), "IntegrinECM_10_sample_example.xlsx")

E <- read_expression_matrix(expression_file, gene_id = "gene")
metadata <- read_sample_metadata(
  sample_metadata_file,
  sample_id = "sample_id",
  group = "analysis_group"
)

result <- score_integrin_ecm(E, metadata)
write_score_workbook(result, output_file)

print(result)
cat("Workbook:", normalizePath(output_file, mustWork = TRUE), "\n")
