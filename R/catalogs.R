#' Parse a stoichiometric structure formula
#'
#' Parses formulas such as `(COL1A1)2(COL1A2)1` into component genes and
#' coefficients.
#'
#' @param formula One character string using repeated `(GENE)coefficient`
#'   terms.
#' @return A data frame with columns `gene` and `coefficient`.
parse_structure_formula <- function(formula) {
  if (!is.character(formula) || length(formula) != 1L || is.na(formula) || formula == "") {
    stop("formula must be one non-empty character string.", call. = FALSE)
  }
  match_position <- gregexpr("\\(([A-Za-z0-9._-]+)\\)([0-9]+)", formula, perl = TRUE)[[1L]]
  if (identical(match_position, -1L)) {
    stop("Could not parse formula: ", formula, call. = FALSE)
  }
  pieces <- regmatches(formula, list(match_position))[[1L]]
  if (!identical(paste0(pieces, collapse = ""), formula)) {
    stop("Formula was only partially parsed: ", formula, call. = FALSE)
  }
  genes <- sub("^\\(([^)]+)\\)[0-9]+$", "\\1", pieces, perl = TRUE)
  coefficients <- as.numeric(sub("^\\([^)]+\\)([0-9]+)$", "\\1", pieces, perl = TRUE))
  if (any(coefficients <= 0)) {
    stop("All formula coefficients must be positive.", call. = FALSE)
  }
  data.frame(gene = genes, coefficient = coefficients, stringsAsFactors = FALSE)
}

.read_catalog_csv <- function(path) {
  if (!file.exists(path)) stop("Catalog file does not exist: ", path, call. = FALSE)
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
}

.catalog_path <- function(filename) {
  path <- system.file("extdata", filename, package = "IntegrinECM")
  if (!nzchar(path)) {
    stop("Bundled catalog is unavailable: ", filename, call. = FALSE)
  }
  path
}

#' Load the curated Integrin-ECM interaction reference
#'
#' @param path Optional path to a replacement interaction-reference CSV.
#' @return A data frame containing the curated interaction records.
load_interaction_reference <- function(path = NULL) {
  if (is.null(path)) path <- .catalog_path("integrin_ecm_reference_284.csv")
  .read_catalog_csv(path)
}

.components_from_formulas <- function(catalog) {
  components <- vector("list", nrow(catalog))
  for (k in seq_len(nrow(catalog))) {
    parsed <- parse_structure_formula(catalog$formula[[k]])
    components[[k]] <- data.frame(
      structure_id = catalog$structure_id[[k]],
      structure_label = catalog$structure_label[[k]],
      gene = parsed$gene,
      coefficient = parsed$coefficient,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, components)
}

.integrin_definition <- function(path) {
  raw <- .read_catalog_csv(path)
  required <- c(
    "integrin", "subunit_formula", "alpha_gene", "beta_gene",
    "alpha_coefficient", "beta_coefficient"
  )
  if (!all(required %in% names(raw))) stop("Invalid Integrin catalog schema.", call. = FALSE)
  catalog <- data.frame(
    structure_id = sprintf("INT%02d", seq_len(nrow(raw))),
    structure_label = raw$integrin,
    formula = raw$subunit_formula,
    alpha_gene = raw$alpha_gene,
    beta_gene = raw$beta_gene,
    alpha_coefficient = as.numeric(raw$alpha_coefficient),
    beta_coefficient = as.numeric(raw$beta_coefficient),
    source_url = if ("source_url" %in% names(raw)) raw$source_url else NA_character_,
    stringsAsFactors = FALSE
  )
  components <- do.call(rbind, lapply(seq_len(nrow(catalog)), function(k) {
    data.frame(
      structure_id = catalog$structure_id[[k]],
      structure_label = catalog$structure_label[[k]],
      gene = c(catalog$alpha_gene[[k]], catalog$beta_gene[[k]]),
      coefficient = c(catalog$alpha_coefficient[[k]], catalog$beta_coefficient[[k]]),
      stringsAsFactors = FALSE
    )
  }))
  list(object = "integrin", catalog = catalog, components = components)
}

.collagen_definition <- function(path) {
  raw <- .read_catalog_csv(path)
  required <- c("catalog_id", "collagen_type", "molecular_species", "formula", "status", "primary_model")
  if (!all(required %in% names(raw))) stop("Invalid Collagen catalog schema.", call. = FALSE)
  primary_model <- as.logical(raw$primary_model)
  keep <- !is.na(primary_model) & primary_model & raw$status == "Included"
  raw <- raw[keep, , drop = FALSE]
  catalog <- data.frame(
    structure_id = sprintf("COL%02d", seq_len(nrow(raw))),
    structure_label = paste0("Collagen ", raw$collagen_type, " [", raw$formula, "]"),
    formula = raw$formula,
    catalog_id = raw$catalog_id,
    collagen_type = raw$collagen_type,
    molecular_species = raw$molecular_species,
    source_url = if ("source_url" %in% names(raw)) raw$source_url else NA_character_,
    evidence_level = if ("evidence_level" %in% names(raw)) raw$evidence_level else NA_character_,
    note = if ("note" %in% names(raw)) raw$note else NA_character_,
    stringsAsFactors = FALSE
  )
  list(object = "collagen", catalog = catalog, components = .components_from_formulas(catalog))
}

.ecm_definition <- function(path, object) {
  raw <- .read_catalog_csv(path)
  required <- c("structure_id", "structure_label", "molecule_class", "subunit_formula")
  if (!all(required %in% names(raw))) stop("Invalid ECM catalog schema.", call. = FALSE)
  source_class <- if (object == "laminin") "Laminin" else "Other fixed ECM"
  raw <- raw[raw$molecule_class == source_class, , drop = FALSE]
  catalog <- data.frame(
    structure_id = raw$structure_id,
    structure_label = raw$structure_label,
    formula = raw$subunit_formula,
    ecm_complex = if ("ecm_complex" %in% names(raw)) raw$ecm_complex else NA_character_,
    assembly_type = if ("assembly_type" %in% names(raw)) raw$assembly_type else NA_character_,
    structure_status = if ("structure_status" %in% names(raw)) raw$structure_status else NA_character_,
    structure_source = if ("structure_source" %in% names(raw)) raw$structure_source else NA_character_,
    structure_url = if ("structure_url" %in% names(raw)) raw$structure_url else NA_character_,
    structure_note = if ("structure_note" %in% names(raw)) raw$structure_note else NA_character_,
    stringsAsFactors = FALSE
  )
  list(object = object, catalog = catalog, components = .components_from_formulas(catalog))
}

#' Load a structure catalog
#'
#' Loads the bundled Integrin, Collagen, Laminin, or Other-ECM structure
#' definition. A replacement raw catalog path may be supplied.
#'
#' @param object One of `"integrin"`, `"collagen"`, `"laminin"`, or
#'   `"other_ecm"`.
#' @param path Optional path to a replacement raw catalog with the corresponding
#'   bundled schema.
#' @return A list containing `object`, structure `catalog`, and long-format
#'   `components`.
load_structure_catalog <- function(
    object = c("integrin", "collagen", "laminin", "other_ecm"),
    path = NULL) {
  object <- match.arg(object)
  if (is.null(path)) {
    path <- switch(
      object,
      integrin = .catalog_path("integrin_heterodimer_catalog.csv"),
      collagen = .catalog_path("human_collagen_structure_catalog.csv"),
      laminin = .catalog_path("ecm_structure_catalog.csv"),
      other_ecm = .catalog_path("ecm_structure_catalog.csv")
    )
  }
  switch(
    object,
    integrin = .integrin_definition(path),
    collagen = .collagen_definition(path),
    laminin = .ecm_definition(path, "laminin"),
    other_ecm = .ecm_definition(path, "other_ecm")
  )
}

.validate_definition <- function(definition) {
  if (!is.list(definition) || !all(c("catalog", "components") %in% names(definition))) {
    stop("definition must contain catalog and components data frames.", call. = FALSE)
  }
  if (!is.data.frame(definition$catalog) || !is.data.frame(definition$components)) {
    stop("definition$catalog and definition$components must be data frames.", call. = FALSE)
  }
  required_catalog <- c("structure_id", "structure_label")
  required_components <- c("structure_id", "structure_label", "gene", "coefficient")
  if (!all(required_catalog %in% names(definition$catalog)) ||
      !all(required_components %in% names(definition$components))) {
    stop("definition has an invalid schema.", call. = FALSE)
  }
  if (anyDuplicated(definition$catalog$structure_id)) {
    stop("Catalog structure_id values must be unique.", call. = FALSE)
  }
  if (anyNA(definition$components$coefficient) || any(definition$components$coefficient <= 0)) {
    stop("Component coefficients must be positive.", call. = FALSE)
  }
  invisible(definition)
}

#' Build a composition matrix
#'
#' @param components Long-format component data frame with `structure_id`,
#'   `gene`, and positive `coefficient` columns.
#' @return Numeric composition matrix with genes in rows and structures in
#'   columns. For NNLS this is manuscript matrix \\eqn{A}; for geometric-mean
#'   scoring its entries are \\eqn{v_{i,k}}.
build_composition_matrix <- function(components) {
  if (!is.data.frame(components) ||
      !all(c("structure_id", "gene", "coefficient") %in% names(components))) {
    stop("components must contain structure_id, gene, and coefficient.", call. = FALSE)
  }
  if (nrow(components) == 0L) stop("components must not be empty.", call. = FALSE)
  if (anyNA(components$coefficient) || any(components$coefficient <= 0)) {
    stop("Component coefficients must be positive.", call. = FALSE)
  }
  if (anyDuplicated(components[c("structure_id", "gene")])) {
    stop("Each structure_id/gene pair must occur once.", call. = FALSE)
  }

  genes <- unique(as.character(components$gene))
  structures <- unique(as.character(components$structure_id))
  A <- matrix(0, length(genes), length(structures), dimnames = list(genes, structures))
  for (row in seq_len(nrow(components))) {
    A[components$gene[[row]], components$structure_id[[row]]] <-
      as.numeric(components$coefficient[[row]])
  }
  A
}
