.validate_component_coverage <- function(E, definition, missing_component_policy, object) {
  missing_component_policy <- .validate_choice(
    missing_component_policy,
    c("error", "omit_structure"),
    "missing_component_policy"
  )
  components <- definition$components
  structure_ids <- unique(as.character(components$structure_id))
  coverage <- lapply(structure_ids, function(k) {
    genes_k <- unique(as.character(components$gene[components$structure_id == k]))
    missing_k <- setdiff(genes_k, rownames(E))
    data.frame(
      structure_id = k,
      structure_label = as.character(components$structure_label[match(k, components$structure_id)]),
      n_required_genes = length(genes_k),
      n_missing = length(missing_k),
      missing_genes = paste(missing_k, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  coverage <- do.call(rbind, coverage)

  missing_rows <- coverage$n_missing > 0
  if (any(missing_rows) && missing_component_policy == "error") {
    stop(
      object, " has structures with missing genes: ",
      paste(coverage$structure_label[missing_rows], collapse = "; "),
      call. = FALSE
    )
  }
  retained <- if (missing_component_policy == "omit_structure") {
    coverage$structure_id[!missing_rows]
  } else {
    coverage$structure_id
  }
  if (length(retained) == 0L) stop("No structures remain after coverage filtering.", call. = FALSE)
  coverage$included <- coverage$structure_id %in% retained

  list(
    catalog = definition$catalog[definition$catalog$structure_id %in% retained, , drop = FALSE],
    components = components[components$structure_id %in% retained, , drop = FALSE],
    audit = coverage
  )
}

.score_matrices_to_long <- function(H, normalized, catalog, metadata, method, object) {
  structure_ids <- rownames(H)
  sample_ids <- colnames(H)
  rows <- vector("list", length(sample_ids))
  labels <- catalog$structure_label[match(structure_ids, catalog$structure_id)]
  groups <- metadata$analysis_group[match(sample_ids, metadata$sample_id)]

  for (s in seq_along(sample_ids)) {
    rows[[s]] <- data.frame(
      method = method,
      object = object,
      sample_id = sample_ids[[s]],
      analysis_group = groups[[s]],
      structure_id = structure_ids,
      structure_label = labels,
      raw_score = as.numeric(H[, s]),
      normalized_fraction = as.numeric(normalized$p[, s]),
      percent_score = as.numeric(normalized$z[, s]),
      rank = as.numeric(normalized$rank[, s]),
      normalization_denominator = as.numeric(normalized$denominator[[s]]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.nnls_fitted_long <- function(E, calculation, metadata) {
  sample_ids <- colnames(E)
  rows <- vector("list", length(sample_ids))
  for (s in seq_along(sample_ids)) {
    rows[[s]] <- data.frame(
      sample_id = sample_ids[[s]],
      analysis_group = metadata$analysis_group[match(sample_ids[[s]], metadata$sample_id)],
      gene = rownames(E),
      observed_expression = as.numeric(E[, s]),
      fitted_expression = as.numeric(calculation$fitted_E[, s]),
      residual = as.numeric(calculation$epsilon[, s]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.geomean_diagnostics <- function(scores) {
  split_scores <- split(scores, scores$sample_id)
  rows <- lapply(split_scores, function(x) {
    denominator <- x$normalization_denominator[[1L]]
    data.frame(
      sample_id = x$sample_id[[1L]],
      analysis_group = x$analysis_group[[1L]],
      n_structures = nrow(x),
      n_zero_scores = sum(x$raw_score == 0, na.rm = TRUE),
      normalization_denominator = denominator,
      normalized_fraction_sum = if (denominator > 0) sum(x$normalized_fraction) else NA_real_,
      percent_score_sum = if (denominator > 0) sum(x$percent_score) else NA_real_,
      all_scores_zero = denominator == 0,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Summarize sample-level scores within groups
#'
#' Group summaries are calculated after sample-level scoring.
#'
#' @param scores Long-format score data returned by `score_structures()`.
#' @param summary_function `"median"` or `"mean"`.
#' @param ties_method Tie method passed to `rank()`.
#' @return A group-level score and rank data frame.
summarize_structure_scores <- function(
    scores,
    summary_function = "median",
    ties_method = "average") {
  summary_function <- .validate_choice(summary_function, c("median", "mean"), "summary_function")
  ties_method <- .validate_choice(
    ties_method,
    c("average", "first", "last", "random", "max", "min"),
    "ties_method"
  )
  required <- c("object", "method", "analysis_group", "structure_id", "structure_label", "percent_score")
  if (!is.data.frame(scores) || !all(required %in% names(scores))) {
    stop("scores has an invalid schema.", call. = FALSE)
  }
  summary_fun <- if (summary_function == "median") stats::median else base::mean
  group_key <- interaction(
    scores$object, scores$method, scores$analysis_group, scores$structure_id,
    drop = TRUE, lex.order = TRUE
  )
  rows <- lapply(split(scores, group_key), function(x) {
    valid <- is.finite(x$percent_score)
    group_score <- if (any(valid)) summary_fun(x$percent_score[valid]) else NA_real_
    data.frame(
      object = x$object[[1L]],
      method = x$method[[1L]],
      analysis_group = x$analysis_group[[1L]],
      structure_id = x$structure_id[[1L]],
      structure_label = x$structure_label[[1L]],
      n_samples = sum(valid),
      summary_function = summary_function,
      group_score = group_score,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rank_key <- interaction(result$object, result$method, result$analysis_group, drop = TRUE)
  result$group_rank <- stats::ave(
    result$group_score,
    rank_key,
    FUN = function(x) rank(-x, ties.method = ties_method, na.last = "keep")
  )
  rownames(result) <- NULL
  result
}

#' Score a predefined set of molecular structures
#'
#' This is the main method-object function. It calculates every sample
#' independently, normalizes scores within the selected structure set, and
#' returns compact scores plus technical diagnostics.
#'
#' @param expression Non-negative numeric matrix with genes in rows and samples
#'   in columns.
#' @param definition Structure definition returned by `load_structure_catalog()`
#'   or a compatible list with `catalog` and `components`.
#' @param method `"nnls"` or `"geomean"`.
#' @param metadata Optional data frame containing `sample_id` and optionally
#'   `analysis_group`.
#' @param object Object label stored in the result.
#' @param missing_component_policy `"error"` or `"omit_structure"`.
#' @param nnls_tolerance Numerical tolerance for NNLS boundary scores.
#' @param output_digits Decimal places for displayed percentage scores.
#' @param rank_ties_method Tie method for within-sample ranks.
#' @param include_zero_scores_in_ranking Whether zero scores receive ranks.
#' @param group_summary_function `"median"` or `"mean"`.
#' @return An `integrin_ecm_module_result`.
score_structures <- function(
    expression,
    definition,
    method = c("nnls", "geomean"),
    metadata = NULL,
    object = NULL,
    missing_component_policy = "error",
    nnls_tolerance = 1e-10,
    output_digits = 2L,
    rank_ties_method = "average",
    include_zero_scores_in_ranking = TRUE,
    group_summary_function = "median") {
  method <- match.arg(method)
  .validate_numeric_matrix(expression, "expression")
  if (is.null(rownames(expression))) stop("expression must have gene row names.", call. = FALSE)
  if (is.null(colnames(expression))) colnames(expression) <- .sample_ids(expression)
  .validate_definition(definition)
  if (is.null(object)) {
    object <- if (!is.null(definition$object)) definition$object else "custom"
  }
  metadata <- .align_metadata(expression, metadata)

  coverage <- .validate_component_coverage(
    expression, definition, missing_component_policy, object
  )
  catalog <- coverage$catalog
  components <- coverage$components
  composition <- build_composition_matrix(components)
  expression_used <- expression[rownames(composition), metadata$sample_id, drop = FALSE]

  if (method == "nnls") {
    calculation <- score_nnls(composition, expression_used, tolerance = nnls_tolerance)
    H <- calculation$x_hat
    method_label <- "NNLS"
  } else {
    calculation <- score_geomean(expression_used, composition)
    H <- calculation$q
    method_label <- "Stoichiometry-weighted geometric mean"
  }

  normalized <- normalize_score_matrix(
    H,
    digits = output_digits,
    ties_method = rank_ties_method,
    include_zero = include_zero_scores_in_ranking
  )
  scores <- .score_matrices_to_long(
    H, normalized, catalog, metadata, method_label, object
  )

  if (method == "nnls") {
    diagnostics <- merge(
      calculation$diagnostics,
      metadata[c("sample_id", "analysis_group")],
      by = "sample_id",
      all.x = TRUE,
      sort = FALSE
    )
    fitted <- .nnls_fitted_long(expression_used, calculation, metadata)
  } else {
    diagnostics <- .geomean_diagnostics(scores)
    fitted <- NULL
  }

  group_summary <- summarize_structure_scores(
    scores,
    summary_function = group_summary_function,
    ties_method = rank_ties_method
  )

  structure(
    list(
      object = object,
      method = method,
      scores = scores,
      group_summary = group_summary,
      diagnostics = diagnostics,
      fitted = fitted,
      catalog = catalog,
      components = components,
      coverage = coverage$audit,
      composition = composition,
      calculation = calculation,
      parameters = list(
        missing_component_policy = missing_component_policy,
        nnls_tolerance = nnls_tolerance,
        output_digits = output_digits,
        rank_ties_method = rank_ties_method,
        include_zero_scores_in_ranking = include_zero_scores_in_ranking,
        group_summary_function = group_summary_function
      ),
      symbol_map = list(
        composition_matrix = if (method == "nnls") "A" else "v_i_k",
        expression = if (method == "nnls") "E_s" else "t_s_i",
        raw_score = if (method == "nnls") "x_hat_s_k" else "q_s_k",
        normalized_fraction = "p_s_k",
        percentage_score = "z_s_k"
      )
    ),
    class = "integrin_ecm_module_result"
  )
}

print.integrin_ecm_module_result <- function(x, ...) {
  cat("IntegrinECM module result\n")
  cat("  object: ", x$object, "\n", sep = "")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  samples: ", length(unique(x$scores$sample_id)), "\n", sep = "")
  cat("  structures: ", length(unique(x$scores$structure_id)), "\n", sep = "")
  invisible(x)
}
