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

.safe_sheet_name <- function(x) {
  substr(gsub("[^A-Za-z0-9_]", "_", x), 1L, 31L)
}

.run_info <- function(result) {
  parameter_values <- vapply(result$parameters, function(x) {
    if (is.null(x)) "NULL" else paste(x, collapse = ";")
  }, character(1L))
  rbind(
    data.frame(
      category = "run",
      item = c("completed_at", "package_version"),
      value = c(
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
}

.set_xlsx_font <- function(path, font_name = "Arial") {
  unpack_dir <- tempfile("integrin_ecm_xlsx_")
  rebuilt_path <- tempfile(fileext = ".xlsx")
  dir.create(unpack_dir, recursive = TRUE)
  on.exit(unlink(c(unpack_dir, rebuilt_path), recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(path, exdir = unpack_dir)
  styles_path <- file.path(unpack_dir, "xl", "styles.xml")
  if (!file.exists(styles_path)) return(invisible(path))
  styles_xml <- paste(readLines(styles_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  styles_xml <- gsub(
    '(<name val=")[^"]+("/>)',
    paste0("\\1", font_name, "\\2"),
    styles_xml,
    perl = TRUE
  )
  styles_xml <- gsub('<scheme val="[^"]+"/>', "", styles_xml, perl = TRUE)
  writeLines(styles_xml, styles_path, useBytes = TRUE)

  files <- list.files(unpack_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  old_directory <- getwd()
  on.exit(setwd(old_directory), add = TRUE)
  setwd(unpack_dir)
  utils::zip(rebuilt_path, files, flags = "-q")
  setwd(old_directory)
  if (!file.copy(rebuilt_path, path, overwrite = TRUE)) {
    stop("Failed to apply workbook font.", call. = FALSE)
  }
  invisible(path)
}

#' Write an Integrin-ECM result workbook
#'
#' The default workbook contains sample metadata, run information, and one
#' compact percentage-score sheet for each selected method-object module.
#'
#' @param result An `integrin_ecm_result`.
#' @param path Output `.xlsx` path.
#' @param optional_sheets Any combination of `"raw_long"`, `"summary"`,
#'   `"diagnostics"`, `"fitted"`, or `"annotations"`.
#' @param font_name Workbook font name.
#' @return Invisibly returns the written sheet names.
write_score_workbook <- function(
    result,
    path,
    optional_sheets = character(0),
    font_name = "Arial") {
  if (!inherits(result, "integrin_ecm_result")) {
    stop("result must be an integrin_ecm_result.", call. = FALSE)
  }
  allowed_optional <- c("raw_long", "summary", "diagnostics", "fitted", "annotations")
  if (any(!optional_sheets %in% allowed_optional)) {
    stop("Unsupported optional_sheets value.", call. = FALSE)
  }
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Package 'writexl' is required to write workbooks.", call. = FALSE)
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
    if ("annotations" %in% optional_sheets) {
      annotation_name <- .safe_sheet_name(paste0(module$object, "_annotations"))
      if (!annotation_name %in% names(sheets)) sheets[[annotation_name]] <- module$catalog
    }
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(sheets, path)
  .set_xlsx_font(path, font_name = font_name)
  invisible(names(sheets))
}
