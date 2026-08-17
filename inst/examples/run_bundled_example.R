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

result <- score_integrin_ecm_files(
  expression_file = expression_file,
  sample_metadata_file = sample_metadata_file,
  gene_id = "gene",
  sample_id = "sample_id",
  group = "analysis_group"
)
write_score_workbook(result, output_file)

print(result)
print(result$input_audit)
cat("Workbook:", normalizePath(output_file, mustWork = TRUE), "\n")
