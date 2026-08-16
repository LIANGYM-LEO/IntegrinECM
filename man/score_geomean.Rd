\name{score_geomean}
\alias{score_geomean}
\alias{score_geomean_structure}
\title{Stoichiometry-weighted geometric-mean scoring}
\usage{
score_geomean_structure(t_s_i, v_i_k)

score_geomean(t, V)
}
\arguments{
  \item{t_s_i}{Non-negative component-expression vector for structure \eqn{k} in sample \eqn{s}.}
  \item{v_i_k}{Stoichiometric coefficient vector.}
  \item{t}{Non-negative gene-by-sample expression matrix.}
  \item{V}{Gene-by-structure stoichiometric matrix.}
}
\description{
Calculates the product \eqn{u_{s,k}} and its total-stoichiometry root
\eqn{q_{s,k}}. No pseudocount is used; a zero required component gives a zero
structure score.
}
\value{
The structure function returns \code{u_s_k} and \code{q_s_k}. The matrix
function returns structures-by-samples matrices \code{u} and \code{q}.
}
