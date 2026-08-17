.scores_to_wide <- function(scores) {
  sample_ids <- unique(scores$sample_id)
  structure_labels <- unique(scores$structure_label)
  z <- matrix(
    NA_real_, length(sample_ids), length(structure_labels),
    dimnames = list(sample_ids, structure_labels)
  )
  for (row in seq_len(nrow(scores))) {
    z[scores$sample_id[[row]], scores$structure_label[[row]]] <- scores$percent_score[[row]]
  }
  groups <- scores$analysis_group[match(sample_ids, scores$sample_id)]
  data.frame(
    sample_id = sample_ids,
    analysis_group = groups,
    z,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

.matrix_to_sample_wide <- function(x) {
  if (!is.matrix(x)) stop("x must be a matrix.", call. = FALSE)
  sample_ids <- colnames(x)
  if (is.null(sample_ids)) sample_ids <- paste0("sample_", seq_len(ncol(x)))
  data.frame(
    sample_id = sample_ids,
    as.data.frame(t(x), check.names = FALSE, stringsAsFactors = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

.composition_to_wide <- function(x) {
  if (!is.matrix(x)) stop("x must be a matrix.", call. = FALSE)
  gene_ids <- rownames(x)
  if (is.null(gene_ids)) gene_ids <- paste0("gene_", seq_len(nrow(x)))
  data.frame(
    gene = gene_ids,
    as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

.safe_sheet_name <- function(x) {
  substr(gsub("[^A-Za-z0-9_]", "_", x), 1L, 31L)
}

.run_info <- function(result) {
  parameter_values <- vapply(result$parameters, function(x) {
    if (is.null(x)) "NULL" else paste(x, collapse = ";")
  }, character(1L))
  run_rows <- rbind(
    data.frame(
      category = "run",
      item = c("run_status", "completed_at", "package_version"),
      value = c(
        if (is.null(result$run_status)) "unknown" else result$run_status,
        format(result$completed_at, "%Y-%m-%d %H:%M:%S %Z"),
        as.character(utils::packageVersion("IntegrinECM"))
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      category = "parameter",
      item = names(parameter_values),
      value = unname(parameter_values),
      stringsAsFactors = FALSE
    ),
    data.frame(
      category = "module",
      item = names(result$modules),
      value = vapply(result$modules, paste, collapse = ";", FUN.VALUE = character(1L)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      category = "output",
      item = "score_definition",
      value = "z_s_k = 100 p_s_k; balanced to sum to 100.00 within each sample and method-object result",
      stringsAsFactors = FALSE
    )
  )

  input_rows <- if (is.null(result$input_audit) || nrow(result$input_audit) == 0L) {
    data.frame(category = character(), item = character(), value = character())
  } else {
    rows <- lapply(seq_len(nrow(result$input_audit)), function(index) {
      audit <- result$input_audit[index, , drop = FALSE]
      prefix <- paste(audit$source_result, audit$input_type, sep = "_")
      path_value <- if (is.na(audit$path) || audit$path == "") audit$input_mode else audit$path
      md5_value <- if (is.na(audit$md5) || audit$md5 == "") "not_available" else audit$md5
      data.frame(
        category = "input",
        item = c(paste0(prefix, "_path"), paste0(prefix, "_md5")),
        value = c(path_value, md5_value),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }

  coverage_rows <- lapply(names(result$results), function(name) {
    coverage <- result$results[[name]]$coverage
    data.frame(
      category = "coverage",
      item = paste0(name, "_structures"),
      value = paste0(
        "included=", sum(coverage$included),
        "; total=", nrow(coverage),
        "; structures_with_missing_genes=", sum(coverage$n_missing > 0)
      ),
      stringsAsFactors = FALSE
    )
  })
  coverage_rows <- do.call(rbind, coverage_rows)

  source_rows <- if (is.null(result$source_runs) || nrow(result$source_runs) == 0L) {
    data.frame(category = character(), item = character(), value = character())
  } else {
    data.frame(
      category = "source_run",
      item = result$source_runs$source_result,
      value = paste0(
        "status=", result$source_runs$run_status,
        "; completed_at=", result$source_runs$completed_at,
        "; modules=", result$source_runs$modules
      ),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, list(run_rows, input_rows, coverage_rows, source_rows))
}

.expand_optional_sheets <- function(optional_sheets) {
  detail_options <- c(
    "raw_long", "summary", "diagnostics", "fitted", "coverage",
    "composition", "geomean_intermediates", "audit", "annotations",
    "components", "interaction_reference"
  )
  aliases <- c("intermediates", "references", "all")
  if (any(!optional_sheets %in% c(detail_options, aliases))) {
    stop("Unsupported optional_sheets value.", call. = FALSE)
  }
  expanded <- optional_sheets
  if (any(optional_sheets %in% c("intermediates", "all"))) {
    expanded <- c(
      expanded,
      "raw_long", "summary", "diagnostics", "fitted", "coverage",
      "composition", "geomean_intermediates", "audit"
    )
  }
  if (any(optional_sheets %in% c("references", "all"))) {
    expanded <- c(expanded, "annotations", "components", "interaction_reference")
  }
  unique(expanded[expanded %in% detail_options])
}

.write_xlsx_safely <- function(sheets, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(sheets, path)
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop("Workbook was not created successfully.", call. = FALSE)
  }
  invisible(names(sheets))
}

.bundled_reference_sheets <- function() {
  sheets <- list(reference_integrin_ecm = load_interaction_reference())
  for (object in c("integrin", "collagen", "laminin", "other_ecm")) {
    definition <- load_structure_catalog(object)
    sheets[[.safe_sheet_name(paste0("ref_", object, "_catalog"))]] <- definition$catalog
    sheets[[.safe_sheet_name(paste0("ref_", object, "_components"))]] <- definition$components
  }
  sheets
}

#' Export all bundled reference resources
#'
#' Writes the curated 284-record Integrin-ECM interaction resource together
#' with all fixed structure catalogs and their component stoichiometry.
#'
#' @param path Output `.xlsx` path.
#' @return Invisibly returns the written sheet names.
write_reference_workbook <- function(path) {
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Package 'writexl' is required to write workbooks.", call. = FALSE)
  }
  .write_xlsx_safely(.bundled_reference_sheets(), path)
}

#' Write an Integrin-ECM result workbook
#'
#' The default workbook contains sample metadata, run information, and one
#' compact percentage-score sheet for each selected method-object module.
#'
#' @param result An `integrin_ecm_result`.
#' @param path Output `.xlsx` path.
#' @param optional_sheets Individual technical-sheet names or the aliases
#'   `"intermediates"`, `"references"`, or `"all"`. See Details.
#' @param font_name Deprecated compatibility argument. Font post-processing is
#'   no longer applied because unpacking and rebuilding XLSX files was not
#'   reliable across operating systems.
#' @return Invisibly returns the written sheet names.
#' @details `"intermediates"` expands to raw long-format scores, group
#'   summaries, diagnostics, fitted values, coverage, composition matrices,
#'   geometric-mean intermediates, and audit tables. `"references"` exports
#'   the selected structure catalogs and component stoichiometry plus the
#'   284-record interaction resource. `"all"` requests both groups.
write_score_workbook <- function(
    result,
    path,
    optional_sheets = character(0),
    font_name = NULL) {
  if (!inherits(result, "integrin_ecm_result")) {
    stop("result must be an integrin_ecm_result.", call. = FALSE)
  }
  optional_sheets <- .expand_optional_sheets(optional_sheets)
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Package 'writexl' is required to write workbooks.", call. = FALSE)
  }
  if (!is.null(font_name)) {
    warning(
      "font_name is deprecated and is not applied; the workbook uses the portable writexl default style.",
      call. = FALSE
    )
  }

  sheets <- list(
    `01_samples` = result$samples,
    `02_run_info` = .run_info(result)
  )
  for (name in names(result$results)) {
    module <- result$results[[name]]
    sheets[[.safe_sheet_name(name)]] <- .scores_to_wide(module$scores)
  }
  for (name in names(result$results)) {
    module <- result$results[[name]]
    if ("raw_long" %in% optional_sheets) {
      sheets[[.safe_sheet_name(paste0(name, "_raw"))]] <- module$scores
    }
    if ("summary" %in% optional_sheets) {
      sheets[[.safe_sheet_name(paste0(name, "_summary"))]] <- module$group_summary
    }
    if ("diagnostics" %in% optional_sheets) {
      sheets[[.safe_sheet_name(paste0(name, "_diagnostics"))]] <- module$diagnostics
    }
    if ("fitted" %in% optional_sheets && !is.null(module$fitted)) {
      sheets[[.safe_sheet_name(paste0(name, "_fitted"))]] <- module$fitted
    }
    if ("coverage" %in% optional_sheets) {
      sheets[[.safe_sheet_name(paste0(name, "_coverage"))]] <- module$coverage
    }
    if ("composition" %in% optional_sheets) {
      sheets[[.safe_sheet_name(paste0(name, "_composition"))]] <-
        .composition_to_wide(module$composition)
    }
    if ("geomean_intermediates" %in% optional_sheets && module$method == "geomean") {
      sheets[[.safe_sheet_name(paste0(name, "_u"))]] <-
        .matrix_to_sample_wide(module$calculation$u)
      sheets[[.safe_sheet_name(paste0(name, "_q"))]] <-
        .matrix_to_sample_wide(module$calculation$q)
    }
    if ("annotations" %in% optional_sheets) {
      annotation_name <- .safe_sheet_name(paste0(module$object, "_annotations"))
      if (!annotation_name %in% names(sheets)) sheets[[annotation_name]] <- module$catalog
    }
    if ("components" %in% optional_sheets) {
      component_name <- .safe_sheet_name(paste0(module$object, "_components"))
      if (!component_name %in% names(sheets)) sheets[[component_name]] <- module$components
    }
  }
  if ("audit" %in% optional_sheets) {
    sheets[["input_audit"]] <- result$input_audit
    if (!is.null(result$source_runs)) sheets[["source_runs"]] <- result$source_runs
  }
  if ("interaction_reference" %in% optional_sheets) {
    sheets[["reference_integrin_ecm"]] <- load_interaction_reference()
  }

  .write_xlsx_safely(sheets, path)
}
