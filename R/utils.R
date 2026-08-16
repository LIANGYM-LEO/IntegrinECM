.validate_numeric_matrix <- function(x, name, nonnegative = TRUE) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(name, " must be a numeric matrix.", call. = FALSE)
  }
  if (length(x) == 0L || anyNA(x) || any(!is.finite(x))) {
    stop(name, " must contain finite, non-missing values.", call. = FALSE)
  }
  if (nonnegative && any(x < 0)) {
    stop(name, " must contain non-negative values.", call. = FALSE)
  }
  invisible(x)
}

.validate_numeric_vector <- function(x, name, nonnegative = TRUE) {
  if (!is.numeric(x) || is.matrix(x) || length(x) == 0L) {
    stop(name, " must be a non-empty numeric vector.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(name, " must contain finite, non-missing values.", call. = FALSE)
  }
  if (nonnegative && any(x < 0)) {
    stop(name, " must contain non-negative values.", call. = FALSE)
  }
  invisible(x)
}

.validate_choice <- function(x, choices, name) {
  if (length(x) != 1L || is.na(x) || !x %in% choices) {
    stop(name, " must be one of: ", paste(choices, collapse = ", "), ".", call. = FALSE)
  }
  x
}

.rank_scores <- function(z_s_k, ties_method = "average", include_zero = TRUE) {
  ties_method <- .validate_choice(
    ties_method,
    c("average", "first", "last", "random", "max", "min"),
    "ties_method"
  )
  valid <- is.finite(z_s_k) & (include_zero | z_s_k != 0)
  rank_s_k <- rep(NA_real_, length(z_s_k))
  rank_s_k[valid] <- rank(-z_s_k[valid], ties.method = ties_method)
  names(rank_s_k) <- names(z_s_k)
  rank_s_k
}

.balanced_percent <- function(h_s_k, digits = 2L) {
  if (length(digits) != 1L || is.na(digits) || digits < 0 || digits != as.integer(digits)) {
    stop("digits must be one non-negative integer.", call. = FALSE)
  }
  if (anyNA(h_s_k) || any(!is.finite(h_s_k)) || any(h_s_k < 0) || sum(h_s_k) <= 0) {
    z_s_k <- rep(NA_real_, length(h_s_k))
    names(z_s_k) <- names(h_s_k)
    return(z_s_k)
  }

  scale_factor <- 10^as.integer(digits)
  exact_units <- h_s_k / sum(h_s_k) * 100 * scale_factor
  displayed_units <- floor(exact_units + 1e-12)
  remainder <- as.integer(round(100 * scale_factor - sum(displayed_units)))

  if (remainder > 0L) {
    priority <- order(-(exact_units - displayed_units), seq_along(h_s_k))
    displayed_units[priority[seq_len(remainder)]] <-
      displayed_units[priority[seq_len(remainder)]] + 1
  }

  z_s_k <- displayed_units / scale_factor
  names(z_s_k) <- names(h_s_k)
  z_s_k
}

.align_expression_to_genes <- function(E, genes, object_name = "expression matrix") {
  .validate_numeric_matrix(E, object_name)
  if (is.null(rownames(E))) {
    stop(object_name, " must have gene identifiers as row names.", call. = FALSE)
  }
  missing_genes <- setdiff(genes, rownames(E))
  if (length(missing_genes) > 0L) {
    stop(
      object_name, " is missing required genes: ",
      paste(utils::head(missing_genes, 20L), collapse = ", "),
      call. = FALSE
    )
  }
  E[genes, , drop = FALSE]
}

.sample_ids <- function(E) {
  ids <- colnames(E)
  if (is.null(ids)) ids <- paste0("sample_", seq_len(ncol(E)))
  if (anyDuplicated(ids)) stop("Sample identifiers must be unique.", call. = FALSE)
  ids
}

.structure_ids <- function(A) {
  ids <- colnames(A)
  if (is.null(ids)) ids <- paste0("structure_", seq_len(ncol(A)))
  if (anyDuplicated(ids)) stop("Structure identifiers must be unique.", call. = FALSE)
  ids
}
