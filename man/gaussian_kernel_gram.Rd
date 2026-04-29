\name{gaussian_kernel_gram}
\alias{gaussian_kernel_gram}
\title{Gaussian Kernel Gram Matrix}
\description{Computes the Gaussian (RBF) kernel Gram matrix for a data matrix.}
\usage{
gaussian_kernel_gram(X, sigma = NULL, normalize_n = FALSE)
}
\arguments{
  \item{X}{an n-by-p numeric matrix.}
  \item{sigma}{bandwidth parameter. If NULL, uses the median heuristic (median of pairwise distances divided by ncol(X)).}
  \item{normalize_n}{logical. If TRUE, normalize the kernel matrix by n. Default is FALSE.}
}
\value{
  An n-by-n symmetric positive semi-definite kernel matrix with attribute "sigma" storing the bandwidth used.
}
\examples{
X <- matrix(rnorm(50), 10, 5)
K <- gaussian_kernel_gram(X)
}
