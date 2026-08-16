#' Read a gene-by-sample expression matrix
#'
#' @param path CSV file path.
#' @param gene_id Name of the gene identifier column.
#' @param duplicate_rule One of `"error"`, `"sum"`, `"mean"`, or `"max"`.
#' @return A non-negative numeric matrix with genes in rows and samples in
#'   columns.
read_expression_matrix <- function(path, gene_id = "gene", duplicate_rule = "error") {
  duplicate_rule <- .validate_choice(
    duplicate_rule, c("error", "sum", "mean", "max"), "duplicate_rule"
  )
  if (!file.exists(path)) stop("Expression file does not exist: ", path, call. = FALSE)
  raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!gene_id %in% names(raw)) stop("Expression gene column not found: ", gene_id, call. = FALSE)
  sample_columns <- setdiff(names(raw), gene_id)
  if (length(sample_columns) == 0L) stop("Expression file contains no sample columns.", call. = FALSE)
  if (anyDuplicated(sample_columns)) stop("Sample identifiers must be unique.", call. = FALSE)

  genes <- as.character(raw[[gene_id]])
  if (anyNA(genes) || any(genes == "")) stop("Gene identifiers must be non-missing.", call. = FALSE)
  values <- as.matrix(raw[sample_columns])
  storage.mode(values) <- "double"
  if (anyNA(values) || any(!is.finite(values))) {
    stop("Expression values must be finite and numeric.", call. = FALSE)
  }
  if (any(values < 0)) stop("Expression values must be non-negative.", call. = FALSE)

  duplicated_genes <- unique(genes[duplicated(genes)])
  if (length(duplicated_genes) > 0L && duplicate_rule == "error") {
    stop(
      "Duplicated genes found: ",
      paste(utils::head(duplicated_genes, 20L), collapse = ", "),
      call. = FALSE
    )
  }
  if (length(duplicated_genes) > 0L) {
    split_rows <- split(seq_along(genes), genes)
    aggregate_one <- switch(
      duplicate_rule,
      sum = function(index) colSums(values[index, , drop = FALSE]),
      mean = function(index) colMeans(values[index, , drop = FALSE]),
      max = function(index) apply(values[index, , drop = FALSE], 2L, max)
    )
    values <- do.call(rbind, lapply(split_rows, aggregate_one))
    genes <- names(split_rows)
  }
  rownames(values) <- genes
  colnames(values) <- sample_columns
  attr(values, "duplicated_genes") <- duplicated_genes
  values
}

#' Read sample metadata
#'
#' @param path CSV file path.
#' @param sample_id Name of the sample identifier column.
#' @param group Name of the analysis-group column. If absent, all samples are
#'   assigned to `"All samples"`.
#' @return A data frame with canonical `sample_id` and `analysis_group` columns.
read_sample_metadata <- function(path, sample_id = "sample_id", group = "analysis_group") {
  if (!file.exists(path)) stop("Metadata file does not exist: ", path, call. = FALSE)
  metadata <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!sample_id %in% names(metadata)) stop("Metadata sample column not found: ", sample_id, call. = FALSE)
  if (!group %in% names(metadata)) metadata[[group]] <- "All samples"
  names(metadata)[names(metadata) == sample_id] <- "sample_id"
  names(metadata)[names(metadata) == group] <- "analysis_group"
  metadata$sample_id <- as.character(metadata$sample_id)
  metadata$analysis_group <- as.character(metadata$analysis_group)
  if (anyDuplicated(metadata$sample_id)) stop("Metadata sample identifiers must be unique.", call. = FALSE)
  metadata
}

.align_metadata <- function(E, metadata = NULL) {
  sample_ids <- .sample_ids(E)
  if (is.null(metadata)) {
    return(data.frame(
      sample_id = sample_ids,
      analysis_group = "All samples",
      stringsAsFactors = FALSE
    ))
  }
  if (!is.data.frame(metadata) || !"sample_id" %in% names(metadata)) {
    stop("metadata must contain sample_id.", call. = FALSE)
  }
  if (!"analysis_group" %in% names(metadata)) metadata$analysis_group <- "All samples"
  if (anyDuplicated(metadata$sample_id)) stop("metadata sample_id values must be unique.", call. = FALSE)
  missing_samples <- setdiff(sample_ids, metadata$sample_id)
  if (length(missing_samples) > 0L) {
    stop("Missing metadata for samples: ", paste(missing_samples, collapse = ", "), call. = FALSE)
  }
  metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]
}
