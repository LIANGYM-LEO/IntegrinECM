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
      input_audit = data.frame(
        source_result = "result_1",
        input_type = c("expression", "sample_metadata"),
        input_mode = "in_memory",
        path = NA_character_,
        md5 = NA_character_,
        stringsAsFactors = FALSE
      ),
      run_status = "completed",
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
  result <- score_integrin_ecm(expression = expression, metadata = metadata, ...)
  input_paths <- c(
    expression = normalizePath(expression_file, winslash = "/", mustWork = TRUE),
    sample_metadata = if (is.null(sample_metadata_file)) {
      NA_character_
    } else {
      normalizePath(sample_metadata_file, winslash = "/", mustWork = TRUE)
    }
  )
  input_md5 <- c(
    expression = unname(tools::md5sum(expression_file)),
    sample_metadata = if (is.null(sample_metadata_file)) {
      NA_character_
    } else {
      unname(tools::md5sum(sample_metadata_file))
    }
  )
  result$input_audit <- data.frame(
    source_result = "result_1",
    input_type = names(input_paths),
    input_mode = "file",
    path = unname(input_paths),
    md5 = unname(input_md5),
    stringsAsFactors = FALSE
  )
  result
}

#' Combine compatible IntegrinECM results
#'
#' Combines separately executed method-object modules after confirming that
#' sample metadata and scoring parameters are identical. When file checksums
#' are available, the expression inputs must also match. Duplicate module names
#' are rejected so that method labels and normalization scopes cannot be
#' silently overwritten.
#'
#' @param ... Two or more `integrin_ecm_result` objects, or one list containing
#'   such objects.
#' @return A combined `integrin_ecm_result`.
combine_integrin_ecm_results <- function(...) {
  inputs <- list(...)
  if (length(inputs) == 1L && is.list(inputs[[1L]]) &&
      !inherits(inputs[[1L]], "integrin_ecm_result")) {
    inputs <- inputs[[1L]]
  }
  if (length(inputs) < 2L ||
      any(!vapply(inputs, inherits, logical(1L), what = "integrin_ecm_result"))) {
    stop("Provide at least two integrin_ecm_result objects.", call. = FALSE)
  }

  source_names <- names(inputs)
  if (is.null(source_names)) source_names <- rep("", length(inputs))
  missing_names <- is.na(source_names) | source_names == ""
  source_names[missing_names] <- paste0("result_", which(missing_names))
  if (anyDuplicated(source_names)) {
    source_names <- make.unique(source_names, sep = "_")
  }

  normalize_frame <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    rownames(x) <- NULL
    x
  }
  reference_samples <- normalize_frame(inputs[[1L]]$samples)
  reference_parameters <- inputs[[1L]]$parameters
  for (index in seq_along(inputs)[-1L]) {
    if (!isTRUE(all.equal(
      reference_samples,
      normalize_frame(inputs[[index]]$samples),
      check.attributes = FALSE
    ))) {
      stop("All results must contain identical sample metadata in the same order.", call. = FALSE)
    }
    if (!identical(reference_parameters, inputs[[index]]$parameters)) {
      stop("All results must use identical scoring and normalization parameters.", call. = FALSE)
    }
  }
  known_expression_md5 <- unlist(lapply(inputs, function(x) {
    audit <- x$input_audit
    if (is.null(audit) || !all(c("input_type", "md5") %in% names(audit))) {
      return(character())
    }
    md5 <- audit$md5[audit$input_type == "expression"]
    md5[!is.na(md5) & md5 != ""]
  }), use.names = FALSE)
  if (length(unique(known_expression_md5)) > 1L) {
    stop("File-based results must use the same expression input checksum.", call. = FALSE)
  }

  result_names <- unlist(lapply(inputs, function(x) names(x$results)), use.names = FALSE)
  if (anyDuplicated(result_names)) {
    duplicates <- unique(result_names[duplicated(result_names)])
    stop(
      "Duplicate method-object modules cannot be combined: ",
      paste(duplicates, collapse = ", "),
      call. = FALSE
    )
  }
  combined_results <- list()
  for (input in inputs) {
    for (result_name in names(input$results)) {
      combined_results[[result_name]] <- input$results[[result_name]]
    }
  }

  combined_modules <- list()
  for (module in combined_results) {
    combined_modules[[module$object]] <- unique(c(
      combined_modules[[module$object]],
      module$method
    ))
  }

  audit_rows <- lapply(seq_along(inputs), function(index) {
    audit <- inputs[[index]]$input_audit
    if (is.null(audit)) {
      audit <- data.frame(
        input_type = c("expression", "sample_metadata"),
        input_mode = "unknown",
        path = NA_character_,
        md5 = NA_character_,
        stringsAsFactors = FALSE
      )
    }
    audit$source_result <- source_names[[index]]
    audit[c("source_result", "input_type", "input_mode", "path", "md5")]
  })

  structure(
    list(
      results = combined_results,
      samples = reference_samples,
      modules = combined_modules,
      parameters = reference_parameters,
      input_audit = do.call(rbind, audit_rows),
      source_runs = data.frame(
        source_result = source_names,
        run_status = vapply(inputs, function(x) {
          if (is.null(x$run_status)) "unknown" else x$run_status
        }, character(1L)),
        completed_at = vapply(inputs, function(x) {
          format(x$completed_at, "%Y-%m-%d %H:%M:%S %Z")
        }, character(1L)),
        modules = vapply(inputs, function(x) {
          paste(names(x$results), collapse = ";")
        }, character(1L)),
        stringsAsFactors = FALSE
      ),
      run_status = "completed",
      completed_at = Sys.time()
    ),
    class = "integrin_ecm_result"
  )
}

print.integrin_ecm_result <- function(x, ...) {
  cat("IntegrinECM result\n")
  cat("  status: ", if (is.null(x$run_status)) "unknown" else x$run_status, "\n", sep = "")
  cat("  samples: ", nrow(x$samples), "\n", sep = "")
  cat("  modules: ", paste(names(x$results), collapse = ", "), "\n", sep = "")
  invisible(x)
}
