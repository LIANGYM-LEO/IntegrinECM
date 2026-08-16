\name{load_structure_catalog}
\alias{load_structure_catalog}
\alias{load_interaction_reference}
\alias{parse_structure_formula}
\alias{build_composition_matrix}
\title{Structure catalog utilities}
\usage{
load_structure_catalog(
  object = c("integrin", "collagen", "laminin", "other_ecm"),
  path = NULL
)

load_interaction_reference(path = NULL)

parse_structure_formula(formula)

build_composition_matrix(components)
}
\arguments{
  \item{object}{Biological object to load.}
  \item{path}{Optional replacement raw-catalog or interaction-reference path.}
  \item{formula}{Formula such as \code{(COL1A1)2(COL1A2)1}.}
  \item{components}{Long component table with structure, gene, and coefficient columns.}
}
\description{
Loads bundled catalogs and the curated interaction reference, parses
stoichiometric formulas, and constructs the gene-by-structure composition
matrix.
}
\value{
\code{load_structure_catalog} returns a definition list;
\code{load_interaction_reference} returns an interaction data frame;
\code{parse_structure_formula} returns a component data frame; and
\code{build_composition_matrix} returns a numeric matrix.
}
