library(IntegrinECM)

# Method 1: exact non-negative reconstruction.
A <- matrix(
  c(1, 0,
    1, 1,
    0, 1),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
)
x_true <- c(S1 = 2, S2 = 3)
E_s <- as.numeric(A %*% x_true)
names(E_s) <- rownames(A)
nnls_one <- score_nnls_sample(A, E_s)
stopifnot(max(abs(nnls_one$x_hat_s - x_true)) < 1e-8)
stopifnot(max(abs(nnls_one$epsilon_s)) < 1e-8)
stopifnot(all(nnls_one$x_hat_s >= 0))
stopifnot(nnls_one$diagnostics$n_active_structures == sum(nnls_one$x_hat_s > 0))

# Method 2: weighted product/root and explicit zero behavior.
geo <- score_geomean_structure(
  t_s_i = c(G1 = 8, G2 = 2),
  v_i_k = c(G1 = 2, G2 = 1)
)
stopifnot(isTRUE(all.equal(unname(geo[["u_s_k"]]), 128)))
stopifnot(isTRUE(all.equal(unname(geo[["q_s_k"]]), 128^(1 / 3))))
geo_zero <- score_geomean_structure(c(8, 0), c(2, 1))
stopifnot(identical(unname(geo_zero), c(0, 0)))

# Normalization: exact displayed closure and all-zero flag behavior.
normalized <- normalize_structure_scores(c(S1 = 1, S2 = 1, S3 = 1), digits = 2)
stopifnot(isTRUE(all.equal(sum(normalized$p_s_k), 1)))
stopifnot(isTRUE(all.equal(sum(normalized$z_s_k), 100)))
stopifnot(identical(normalized$z_s_k, c(33.34, 33.33, 33.33)))
all_zero <- normalize_structure_scores(c(S1 = 0, S2 = 0))
stopifnot(all(is.na(all_zero$p_s_k)))
stopifnot(all(is.na(all_zero$z_s_k)))
stopifnot(all(is.na(all_zero$rank_s_k)))

# General method-object interface.
catalog <- data.frame(
  structure_id = c("S1", "S2"),
  structure_label = c("Structure 1", "Structure 2"),
  stringsAsFactors = FALSE
)
components <- data.frame(
  structure_id = c("S1", "S1", "S2", "S2"),
  structure_label = c("Structure 1", "Structure 1", "Structure 2", "Structure 2"),
  gene = c("G1", "G2", "G2", "G3"),
  coefficient = c(1, 1, 1, 1),
  stringsAsFactors = FALSE
)
definition <- list(object = "custom", catalog = catalog, components = components)
expression <- matrix(
  c(5, 7,
    8, 6,
    3, 4),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("G1", "G2", "G3"), c("sample_A", "sample_B"))
)
metadata <- data.frame(
  sample_id = c("sample_A", "sample_B"),
  analysis_group = c("A", "B"),
  stringsAsFactors = FALSE
)
module_nnls <- score_structures(expression, definition, "nnls", metadata, object = "custom")
module_geo <- score_structures(expression, definition, "geomean", metadata, object = "custom")
stopifnot(inherits(module_nnls, "integrin_ecm_module_result"))
stopifnot(inherits(module_geo, "integrin_ecm_module_result"))
stopifnot(all(tapply(module_nnls$scores$percent_score, module_nnls$scores$sample_id, sum) == 100))
stopifnot(all(tapply(module_geo$scores$percent_score, module_geo$scores$sample_id, sum) == 100))

# Bundled catalogs remain readable and retain expected core dimensions.
integrin <- load_structure_catalog("integrin")
collagen <- load_structure_catalog("collagen")
interactions <- load_interaction_reference()
stopifnot(nrow(integrin$catalog) == 24L)
stopifnot(nrow(collagen$catalog) == 33L)
stopifnot(nrow(interactions) == 284L)
stopifnot(all(c("structure_id", "gene", "coefficient") %in% names(collagen$components)))
stopifnot(identical(
  collagen$catalog$structure_label[[1L]],
  "Collagen I [(COL1A1)2(COL1A2)1]"
))

# The bundled example inputs are strictly paired and run through all defaults.
example_expression_path <- system.file(
  "extdata", "example_expression_10_samples.csv",
  package = "IntegrinECM"
)
example_metadata_path <- system.file(
  "extdata", "example_sample_metadata_10_samples.csv",
  package = "IntegrinECM"
)
stopifnot(nzchar(example_expression_path), nzchar(example_metadata_path))
example_expression <- read_expression_matrix(example_expression_path)
example_metadata <- read_sample_metadata(example_metadata_path)
stopifnot(ncol(example_expression) == 10L)
stopifnot(nrow(example_metadata) == 10L)
stopifnot(identical(colnames(example_expression), example_metadata$sample_id))

example_result <- score_integrin_ecm_files(
  example_expression_path,
  example_metadata_path
)
stopifnot(identical(
  names(example_result$results),
  c("integrin_nnls", "collagen_nnls", "laminin_geomean", "other_ecm_geomean")
))
stopifnot(identical(example_result$run_status, "completed"))
stopifnot(all(example_result$input_audit$input_mode == "file"))
stopifnot(all(file.exists(example_result$input_audit$path)))
stopifnot(all(nchar(example_result$input_audit$md5) == 32L))
example_workbook <- tempfile(fileext = ".xlsx")
example_sheets <- write_score_workbook(example_result, example_workbook)
stopifnot(file.exists(example_workbook), file.info(example_workbook)$size > 0)
stopifnot(identical(
  example_sheets,
  c("01_samples", "02_run_info", names(example_result$results))
))
unlink(example_workbook)

# Every intermediate and selected reference can be exported on request.
complete_workbook <- tempfile(fileext = ".xlsx")
complete_sheets <- write_score_workbook(
  example_result,
  complete_workbook,
  optional_sheets = "all"
)
stopifnot(file.exists(complete_workbook), file.info(complete_workbook)$size > 0)
stopifnot(all(c(
  "input_audit",
  "reference_integrin_ecm",
  "collagen_nnls_coverage",
  "collagen_nnls_composition",
  "laminin_geomean_u",
  "laminin_geomean_q",
  "collagen_annotations",
  "collagen_components"
) %in% complete_sheets))
unlink(complete_workbook)

# All bundled reference resources can be exported without a scoring run.
reference_workbook <- tempfile(fileext = ".xlsx")
reference_sheets <- write_reference_workbook(reference_workbook)
stopifnot(file.exists(reference_workbook), file.info(reference_workbook)$size > 0)
stopifnot(length(reference_sheets) == 9L)
stopifnot(all(c(
  "reference_integrin_ecm",
  "ref_integrin_catalog",
  "ref_integrin_components",
  "ref_collagen_catalog",
  "ref_collagen_components",
  "ref_laminin_catalog",
  "ref_laminin_components",
  "ref_other_ecm_catalog",
  "ref_other_ecm_components"
) %in% reference_sheets))
unlink(reference_workbook)

# Separately executed compatible modules combine without changing module names.
integrin_only <- score_integrin_ecm(
  example_expression,
  example_metadata,
  modules = list(integrin = "nnls")
)
collagen_only <- score_integrin_ecm(
  example_expression,
  example_metadata,
  modules = list(collagen = "nnls")
)
combined <- combine_integrin_ecm_results(
  integrin_run = integrin_only,
  collagen_run = collagen_only
)
stopifnot(identical(names(combined$results), c("integrin_nnls", "collagen_nnls")))
stopifnot(identical(combined$modules, list(integrin = "nnls", collagen = "nnls")))
stopifnot(nrow(combined$source_runs) == 2L)
stopifnot(all(combined$input_audit$source_result %in% c("integrin_run", "collagen_run")))

duplicate_error <- try(
  combine_integrin_ecm_results(integrin_only, integrin_only),
  silent = TRUE
)
stopifnot(inherits(duplicate_error, "try-error"))

incompatible <- collagen_only
incompatible$parameters$output_digits <- 3L
parameter_error <- try(
  combine_integrin_ecm_results(integrin_only, incompatible),
  silent = TRUE
)
stopifnot(inherits(parameter_error, "try-error"))

incompatible_samples <- collagen_only
incompatible_samples$samples$analysis_group[[1L]] <- "Different group"
sample_error <- try(
  combine_integrin_ecm_results(integrin_only, incompatible_samples),
  silent = TRUE
)
stopifnot(inherits(sample_error, "try-error"))

file_integrin <- example_result
file_integrin$results <- file_integrin$results["integrin_nnls"]
file_collagen <- example_result
file_collagen$results <- file_collagen$results["collagen_nnls"]
file_collagen$input_audit$md5[file_collagen$input_audit$input_type == "expression"] <-
  paste(rep("0", 32L), collapse = "")
checksum_error <- try(
  combine_integrin_ecm_results(file_integrin, file_collagen),
  silent = TRUE
)
stopifnot(inherits(checksum_error, "try-error"))
