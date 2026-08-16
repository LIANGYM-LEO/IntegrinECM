#' Score one sample by non-negative least squares
#'
#' Implements the manuscript model
#' \\eqn{E_s = A x_s + \\varepsilon_s} and estimates
#' \\eqn{\\hat{x}_s = \\arg\\min_{x_s \\ge 0} ||A x_s - E_s||_2^2}.
#'
#' @param A Numeric composition matrix with component genes in rows and
#'   candidate structures in columns.
#' @param E_s Non-negative expression vector for sample \\eqn{s}, ordered as
#'   the rows of `A`. Named values are aligned to `rownames(A)`.
#' @param tolerance Values of \\eqn{\\hat{x}_s} smaller than this threshold are
#'   set to zero.
#' @return An object containing `x_hat_s`, fitted expression, residual
#'   `epsilon_s`, and diagnostics.
score_nnls_sample <- function(A, E_s, tolerance = 1e-10) {
  .validate_numeric_matrix(A, "A")
  .validate_numeric_vector(E_s, "E_s")
  if (nrow(A) != length(E_s)) {
    stop("nrow(A) must equal length(E_s).", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || tolerance < 0) {
    stop("tolerance must be one non-negative number.", call. = FALSE)
  }

  if (!is.null(rownames(A)) && !is.null(names(E_s))) {
    missing_genes <- setdiff(rownames(A), names(E_s))
    if (length(missing_genes) > 0L) {
      stop("E_s is missing names present in rownames(A).", call. = FALSE)
    }
    E_s <- E_s[rownames(A)]
  }

  fit <- nnls::nnls(A, E_s)
  x_hat_s <- as.numeric(stats::coef(fit))
  x_hat_s[abs(x_hat_s) < tolerance] <- 0
  names(x_hat_s) <- .structure_ids(A)

  fitted_E_s <- as.numeric(A %*% x_hat_s)
  names(fitted_E_s) <- rownames(A)
  epsilon_s <- as.numeric(E_s - fitted_E_s)
  names(epsilon_s) <- rownames(A)

  residual_norm <- sqrt(sum(epsilon_s^2))
  E_s_norm <- sqrt(sum(E_s^2))

  structure(
    list(
      x_hat_s = x_hat_s,
      fitted_E_s = fitted_E_s,
      epsilon_s = epsilon_s,
      diagnostics = list(
        matrix_rank = qr(A, tol = tolerance)$rank,
        n_component_genes = nrow(A),
        n_candidate_structures = ncol(A),
        condition_number = tryCatch(kappa(A), error = function(e) NA_real_),
        residual_norm = residual_norm,
        relative_residual_norm = if (E_s_norm == 0) NA_real_ else residual_norm / E_s_norm,
        n_boundary_zero_scores = sum(x_hat_s == 0),
        solver_mode = fit$mode,
        n_active_structures = fit$nsetp
      )
    ),
    class = "integrin_ecm_nnls_sample"
  )
}

#' Score multiple samples by non-negative least squares
#'
#' Applies `score_nnls_sample()` independently to every sample column.
#'
#' @param A Numeric composition matrix \\eqn{A} with genes by structures.
#' @param E Numeric expression matrix with genes by samples.
#' @param tolerance Non-negative numerical tolerance.
#' @return A list containing matrices `x_hat`, `fitted_E`, `epsilon`, and a
#'   sample-level diagnostics data frame.
score_nnls <- function(A, E, tolerance = 1e-10) {
  .validate_numeric_matrix(A, "A")
  .validate_numeric_matrix(E, "E")

  if (!is.null(rownames(A))) {
    E <- .align_expression_to_genes(E, rownames(A), "E")
  } else if (nrow(A) != nrow(E)) {
    stop("A and E must have the same number of rows.", call. = FALSE)
  }

  sample_ids <- .sample_ids(E)
  structure_ids <- .structure_ids(A)
  gene_ids <- rownames(A)
  if (is.null(gene_ids)) gene_ids <- paste0("gene_", seq_len(nrow(A)))

  x_hat <- matrix(
    0, nrow = ncol(A), ncol = ncol(E),
    dimnames = list(structure_ids, sample_ids)
  )
  fitted_E <- matrix(
    0, nrow = nrow(A), ncol = ncol(E),
    dimnames = list(gene_ids, sample_ids)
  )
  epsilon <- fitted_E
  diagnostics <- vector("list", ncol(E))

  for (s in seq_len(ncol(E))) {
    E_s <- E[, s]
    names(E_s) <- gene_ids
    scored_s <- score_nnls_sample(A, E_s, tolerance = tolerance)
    x_hat[, s] <- scored_s$x_hat_s
    fitted_E[, s] <- scored_s$fitted_E_s
    epsilon[, s] <- scored_s$epsilon_s
    diagnostics[[s]] <- data.frame(
      sample_id = sample_ids[[s]],
      matrix_rank = scored_s$diagnostics$matrix_rank,
      n_component_genes = scored_s$diagnostics$n_component_genes,
      n_candidate_structures = scored_s$diagnostics$n_candidate_structures,
      condition_number = scored_s$diagnostics$condition_number,
      residual_norm = scored_s$diagnostics$residual_norm,
      relative_residual_norm = scored_s$diagnostics$relative_residual_norm,
      n_boundary_zero_scores = scored_s$diagnostics$n_boundary_zero_scores,
      solver_mode = scored_s$diagnostics$solver_mode,
      n_active_structures = scored_s$diagnostics$n_active_structures,
      stringsAsFactors = FALSE
    )
  }

  list(
    x_hat = x_hat,
    fitted_E = fitted_E,
    epsilon = epsilon,
    diagnostics = do.call(rbind, diagnostics),
    A = A
  )
}
