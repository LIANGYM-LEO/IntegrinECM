#' Normalize and rank one sample's structure scores
#'
#' For raw structure scores \\eqn{h_{s,k}}, calculates
#' \\eqn{p_{s,k}=h_{s,k}/\\sum_{l\\in S}h_{s,l}} and
#' \\eqn{z_{s,k}=100p_{s,k}}. Displayed percentage scores use balanced
#' largest-remainder rounding so valid scores sum exactly to 100 at the
#' requested precision.
#'
#' @param h_s_k Non-negative vector of raw structure scores for one sample.
#' @param digits Number of decimal places for `z_s_k`.
#' @param ties_method Tie method passed to `rank()`.
#' @param include_zero Whether zero percentage scores receive ranks.
#' @return A data frame containing `structure_id`, `h_s_k`, `p_s_k`, `z_s_k`,
#'   and `rank_s_k`.
normalize_structure_scores <- function(
    h_s_k,
    digits = 2L,
    ties_method = "average",
    include_zero = TRUE) {
  .validate_numeric_vector(h_s_k, "h_s_k")
  if (!is.logical(include_zero) || length(include_zero) != 1L || is.na(include_zero)) {
    stop("include_zero must be TRUE or FALSE.", call. = FALSE)
  }

  structure_id <- names(h_s_k)
  if (is.null(structure_id)) structure_id <- paste0("structure_", seq_along(h_s_k))
  denominator <- sum(h_s_k)

  if (denominator <= 0) {
    p_s_k <- z_s_k <- rank_s_k <- rep(NA_real_, length(h_s_k))
  } else {
    p_s_k <- h_s_k / denominator
    names(p_s_k) <- structure_id
    z_s_k <- .balanced_percent(h_s_k, digits = digits)
    rank_s_k <- .rank_scores(z_s_k, ties_method = ties_method, include_zero = include_zero)
  }

  data.frame(
    structure_id = structure_id,
    h_s_k = as.numeric(h_s_k),
    p_s_k = as.numeric(p_s_k),
    z_s_k = as.numeric(z_s_k),
    rank_s_k = as.numeric(rank_s_k),
    normalization_denominator = denominator,
    stringsAsFactors = FALSE
  )
}

#' Normalize a structures-by-samples score matrix
#'
#' @param H Non-negative numeric matrix with structures in rows and samples in
#'   columns. Column \\eqn{s} is the vector \\eqn{h_{s,k}} over \\eqn{k}.
#' @inheritParams normalize_structure_scores
#' @return A list with matrices `p`, `z`, `rank`, and sample denominators.
normalize_score_matrix <- function(
    H,
    digits = 2L,
    ties_method = "average",
    include_zero = TRUE) {
  .validate_numeric_matrix(H, "H")
  structure_ids <- rownames(H)
  if (is.null(structure_ids)) structure_ids <- paste0("structure_", seq_len(nrow(H)))
  sample_ids <- .sample_ids(H)

  p <- z <- rank_matrix <- matrix(
    NA_real_, nrow(H), ncol(H),
    dimnames = list(structure_ids, sample_ids)
  )
  denominator <- stats::setNames(numeric(ncol(H)), sample_ids)

  for (s in seq_len(ncol(H))) {
    h_s_k <- H[, s]
    names(h_s_k) <- structure_ids
    normalized_s <- normalize_structure_scores(
      h_s_k,
      digits = digits,
      ties_method = ties_method,
      include_zero = include_zero
    )
    p[, s] <- normalized_s$p_s_k
    z[, s] <- normalized_s$z_s_k
    rank_matrix[, s] <- normalized_s$rank_s_k
    denominator[[s]] <- normalized_s$normalization_denominator[[1L]]
  }

  list(p = p, z = z, rank = rank_matrix, denominator = denominator)
}
