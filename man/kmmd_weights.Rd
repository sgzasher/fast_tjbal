\name{kmmd_weights}
\alias{kmmd_weights}
\title{Kernel MMD Balancing Weights}
\description{Computes balancing weights by minimizing the squared kernel Maximum Mean Discrepancy (MMD) between treated and control groups with entropy regularization.}
\usage{
kmmd_weights(K, treated_idx, control_idx, eta = 0.01, max_iter = 5000, tol = 1e-6)
}
\arguments{
  \item{K}{an n-by-n kernel Gram matrix (e.g., from \code{\link{gaussian_kernel_gram}}).}
  \item{treated_idx}{integer vector of row/column indices for treated units.}
  \item{control_idx}{integer vector of row/column indices for control units.}
  \item{eta}{entropy regularization parameter. Larger values produce more uniform weights. Default is 0.01.}
  \item{max_iter}{maximum number of exponentiated gradient iterations. Default is 5000.}
  \item{tol}{convergence tolerance on max absolute change in weights. Default is 1e-6.}
}
\value{
  A list with components:
  \item{w}{numeric vector of weights summing to 1.}
  \item{objective}{the final objective value.}
  \item{converged}{logical indicating whether the solver converged.}
  \item{iterations}{number of iterations performed.}
}
\examples{
X <- matrix(rnorm(100), 20, 5)
K <- gaussian_kernel_gram(X)
sol <- kmmd_weights(K, treated_idx = 1:5, control_idx = 6:20, eta = 0.1)
}
