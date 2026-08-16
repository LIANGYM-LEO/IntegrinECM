\name{normalize_structure_scores}
\alias{normalize_structure_scores}
\alias{normalize_score_matrix}
\title{Within-sample normalization, percentage scoring, and ranking}
\usage{
normalize_structure_scores(
  h_s_k,
  digits = 2L,
  ties_method = "average",
  include_zero = TRUE
)

normalize_score_matrix(
  H,
  digits = 2L,
  ties_method = "average",
  include_zero = TRUE
)
}
\arguments{
  \item{h_s_k}{Non-negative raw structure scores for one sample.}
  \item{H}{Structures-by-samples raw score matrix.}
  \item{digits}{Displayed percentage-score decimal places.}
  \item{ties_method}{Tie method passed to \code{rank}.}
  \item{include_zero}{Whether zero scores receive ranks.}
}
\description{
Computes \eqn{p_{s,k}=h_{s,k}/\sum_{l\in S}h_{s,l}} and
\eqn{z_{s,k}=100p_{s,k}}. Balanced rounding makes each valid displayed sample
total exactly 100 at the requested precision.
}
\value{
The vector function returns a data frame. The matrix function returns matrices
\code{p}, \code{z}, and \code{rank}, plus sample denominators.
}
