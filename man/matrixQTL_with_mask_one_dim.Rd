% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rlib_matrix_ls_with_mask.R
\name{matrixQTL_with_mask_one_dim}
\alias{matrixQTL_with_mask_one_dim}
\title{Solve Y = Xb + e with mask and weight in matrix form}
\usage{
matrixQTL_with_mask_one_dim(Y, X, M, Weight = NULL)
}
\arguments{
\item{Y}{response to regress against (dimension = N x K)}

\item{X}{P predictors to perform regression separately (dimension = N x P)}

\item{M}{mask for regression of Y_k (dimension = N x K)}

\item{Weight}{weights for regression of Y_k (dimension = N x K; default = NULL)}
}
\value{
a list of summary statistics
        beta_hat: estimated b, b_hat (dimension = K x P)
        beta_se: standard deviation of b_hat (dimension = K x P)
}
\description{
For each column i in X, solve Y_k = X_i b_i + e with M_k as mask (entries with M_k = 0 are excluded)
and observations are weighted by W_k
as least squares problem and output estimated effect size and standard deviation
}
\examples{
matrixQTL_with_mask_one_dim(
  Y = matrix(rnorm(300), ncol = 3),
  X = matrix(sample(c(0, 0.5, 1), 200, replace = TRUE), ncol = 2),
  M = matrix(sample(c(0, 1), 300, replace = TRUE), ncol = 1),
  Weight = matrix(runif(300), ncol = 3)
)

}
