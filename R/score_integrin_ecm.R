#' Run modular Integrin-ECM scoring
#'
#' Runs any selected method-object combinations. By default, Integrin and
#' Collagen use NNLS, whereas Laminin and Other ECM use the weighted geometric
#' mean.
#'
#' @param expression Non-negative gene-by-sample numeric matrix.
#' @param metadata Optional sample metadata.
#' @param modules Named list whose entries contain `"nnls"`, `"geomean"`,
#'   both, or `character(0)`.
#' @param catalogs Optional named list of structure definitions. Missing entries
#'   are loaded from the bundled catalogs.
#' @inheritParams score_structures
#' @return An `integrin_ecm_result` containing all selected module results.
score_integrin_ecm <- function(
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
    group_summary_function = "median") {
  .validate_numeric_matrix(expression, "expression")
  if (is.null(rownames(expression))) stop("expression must have gene row names.", call. = FALSE)
  if (is.null(colnames(expression))) colnames(expression) <- .sample_ids(expression)
  metadata <- .align_metadata(expression, metadata)

  supported_objects <- c("integrin", "collagen", "laminin", "other_ecm")
  if (!is.list(modules) || is.null(names(modules)) || any(!names(modules) %in% supported_objects)) {
    stop("modules must be a named list containing supported objects.", call. = FALSE)
  }
  if (is.null(catalogs)) catalogs <- list()
  if (!is.list(catalogs)) stop("catalogs must be NULL or a named list.", call. = FALSE)

  results <- list()
  for (object in names(modules)) {
    methods <- modules[[object]]
    if (length(methods) == 0L) next
    if (any(!methods %in% c("nnls", "geomean"))) {
      stop("Unsupported method for ", object, ".", call. = FALSE)
    }
    if (object %in% c("laminin", "other_ecm") && "nnls" %in% methods) {
      stop("NNLS is not enabled for ", object, "; use geomean.", call. = FALSE)
    }
    definition <- catalogs[[object]]
    if (is.null(definition)) definition <- load_structure_catalog(object)

    for (method in methods) {
      result_name <- paste(object, method, sep = "_")
      results[[result_name]] <- score_structures(
        expression = expression,
        definition = definition,
        method = method,
        metadata = metadata,
        object = object,
        missing_component_policy = missing_component_policy,
        nnls_tolerance = nnls_tolerance,
        output_digits = output_digits,
        rank_ties_method = rank_ties_method,
        include_zero_scores_in_ranking = include_zero_scores_in_ranking,
        group_summary_function = group_summary_function
      )
    }
  }
  if (length(results) == 0L) stop("No method-object modules were selected.", call. = FALSE)

  structure(
    list(
      results = results,
      samples = metadata,
      modules = modules,
      parameters = list(
        missing_component_policy = missing_component_policy,
        nnls_tolerance = nnls_tolerance,
        output_digits = output_digits,
        rank_ties_method = rank_ties_method,
        include_zero_scores_in_ranking = include_zero_scores_in_ranking,
        group_summary_function = group_summary_function
      ),
      completed_at = Sys.time()
    ),
    class = "integrin_ecm_result"
  )
}

#' Run modular scoring from CSV files
#'
#' @param expression_file Gene-by-sample expression CSV.
#' @param sample_metadata_file Optional metadata CSV.
#' @param gene_id Expression gene-column name.
#' @param sample_id Metadata sample-column name.
#' @param group Metadata group-column name.
#' @param duplicate_rule Duplicate-gene rule.
#' @param ... Additional arguments passed to `score_integrin_ecm()`.
#' @return An `integrin_ecm_result`.
score_integrin_ecm_files <- function(
    expression_file,
    sample_metadata_file = NULL,
    gene_id = "gene",
    sample_id = "sample_id",
    group = "analysis_group",
    duplicate_rule = "error",
    ...) {
  expression <- read_expression_matrix(expression_file, gene_id, duplicate_rule)
  metadata <- if (is.null(sample_metadata_file)) NULL else {
    read_sample_metadata(sample_metadata_file, sample_id, group)
  }
  score_integrin_ecm(expression = expression, metadata = metadata, ...)
}

print.integrin_ecm_result <- function(x, ...) {
  cat("IntegrinECM result\n")
  cat("  samples: ", nrow(x$samples), "\n", sep = "")
  cat("  modules: ", paste(names(x$results), collapse = ", "), "\n", sep = "")
  invisible(x)
}
