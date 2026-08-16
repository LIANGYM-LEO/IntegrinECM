#' Score one structure by the stoichiometry-weighted geometric mean
#'
#' Implements \\eqn{u_{s,k}=\\prod_{i\\in I_k}(t_{s,i})^{v_{i,k}}} and
#' \\eqn{q_{s,k}=u_{s,k}^{1/\\sum_{i\\in I_k}v_{i,k}}}. No pseudocount is
#' added. If any required component has zero expression, both scores are zero.
#'
#' @param t_s_i Non-negative expression vector for the component genes of
#'   structure \\eqn{k} in sample \\eqn{s}.
#' @param v_i_k Non-negative stoichiometric coefficient vector aligned with
#'   `t_s_i`. At least one coefficient must be positive.
#' @return A named numeric vector containing `u_s_k` and `q_s_k`.
score_geomean_structure <- function(t_s_i, v_i_k) {
  .validate_numeric_vector(t_s_i, "t_s_i")
  .validate_numeric_vector(v_i_k, "v_i_k")
  if (length(t_s_i) != length(v_i_k)) {
    stop("t_s_i and v_i_k must have the same length.", call. = FALSE)
  }

  if (!is.null(names(t_s_i)) && !is.null(names(v_i_k))) {
    missing_genes <- setdiff(names(v_i_k), names(t_s_i))
    if (length(missing_genes) > 0L) {
      stop("t_s_i is missing genes named in v_i_k.", call. = FALSE)
    }
    t_s_i <- t_s_i[names(v_i_k)]
  }

  required <- v_i_k > 0
  if (!any(required)) {
    stop("v_i_k must contain at least one positive coefficient.", call. = FALSE)
  }
  t_s_i <- t_s_i[required]
  v_i_k <- v_i_k[required]

  if (any(t_s_i == 0)) {
    return(c(u_s_k = 0, q_s_k = 0))
  }

  log_u_s_k <- sum(v_i_k * log(t_s_i))
  total_stoichiometry <- sum(v_i_k)
  u_s_k <- if (log_u_s_k > log(.Machine$double.xmax)) Inf else exp(log_u_s_k)
  q_s_k <- exp(log_u_s_k / total_stoichiometry)
  c(u_s_k = u_s_k, q_s_k = q_s_k)
}

#' Score structures and samples by weighted geometric means
#'
#' @param t Numeric expression matrix with component genes in rows and samples
#'   in columns. The stored element `t[i, s]` represents manuscript quantity
#'   \\eqn{t_{s,i}}.
#' @param V Numeric stoichiometric matrix with component genes in rows and
#'   structures in columns. `V[i, k]` represents \\eqn{v_{i,k}}.
#' @return A list with unrooted product matrix `u` and scale-adjusted score
#'   matrix `q`, both arranged as structures by samples.
score_geomean <- function(t, V) {
  .validate_numeric_matrix(t, "t")
  .validate_numeric_matrix(V, "V")

  if (!is.null(rownames(V))) {
    t <- .align_expression_to_genes(t, rownames(V), "t")
  } else if (nrow(t) != nrow(V)) {
    stop("t and V must have the same number of rows.", call. = FALSE)
  }

  sample_ids <- .sample_ids(t)
  structure_ids <- .structure_ids(V)
  u <- matrix(NA_real_, ncol(V), ncol(t), dimnames = list(structure_ids, sample_ids))
  q <- u

  for (k in seq_len(ncol(V))) {
    v_i_k <- V[, k]
    for (s in seq_len(ncol(t))) {
      t_s_i <- t[, s]
      scored_s_k <- score_geomean_structure(t_s_i, v_i_k)
      u[k, s] <- scored_s_k[["u_s_k"]]
      q[k, s] <- scored_s_k[["q_s_k"]]
    }
  }

  list(u = u, q = q, V = V)
}
