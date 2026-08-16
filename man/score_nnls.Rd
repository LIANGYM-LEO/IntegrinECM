\name{score_nnls}
\alias{score_nnls}
\alias{score_nnls_sample}
\title{Non-negative least-squares structure allocation}
\usage{
score_nnls_sample(A, E_s, tolerance = 1e-10)

score_nnls(A, E, tolerance = 1e-10)
}
\arguments{
  \item{A}{Numeric composition matrix with genes in rows and structures in columns.}
  \item{E_s}{Non-negative expression vector for sample \eqn{s}.}
  \item{E}{Non-negative gene-by-sample expression matrix.}
  \item{tolerance}{Values of the estimated score below this threshold are set to zero.}
}
\description{
Implements \eqn{E_s=A x_s+\varepsilon_s} and estimates the non-negative vector
\eqn{\hat{x}_s} that minimizes \eqn{\lVert A x_s-E_s\rVert_2^2}.
}
\value{
The sample function returns \code{x_hat_s}, fitted expression,
\code{epsilon_s}, and diagnostics. The matrix function returns corresponding
structure-by-sample and gene-by-sample matrices.
}
