\name{tjbal}
\alias{tjbal}
\alias{tjbal.default}
\alias{tjbal.formula}
\title{Trajectory Balancing}
\description{Implements kernel balancing via KMMD minimization with time-series cross-sectional data.}
\usage{
tjbal(\dots)

\method{tjbal}{formula}(formula, data, X.avg.time = NULL, index,
  trim.npre = 0, Y.match.time = NULL, Y.match.npre = NULL,
  demean = TRUE, estimator = "kernel", sigma = NULL, eta = 0.01,
  nystrom = NULL, nystrom_m = 200L, print.baltable = TRUE,
  vce = "jackknife", conf.lvl = 0.95, nsims = 200,
  parallel = TRUE, cores = NULL, seed = NULL, \dots)

\method{tjbal}{default}(data, Y, D, X = NULL, X.avg.time = NULL,
  index, trim.npre = 0, Y.match.time = NULL, Y.match.npre = NULL,
  demean = TRUE, estimator = "kernel", sigma = NULL, eta = 0.01,
  nystrom = NULL, nystrom_m = 200L, print.baltable = TRUE,
  vce = "jackknife", conf.lvl = 0.95, nsims = 200,
  parallel = TRUE, cores = NULL, seed = NULL, \dots)
}
\arguments{
\item{formula}{an object of class "formula": a symbolic description of
  the model to be fitted. The first variable on the right-hand-side is a dichotomous treatment indicator (D); the rest of the right-hand-side variables are time-invariant controls (X). If X are time-varying,
  they will be averaged within each unit based on values specified by \code{X.avg.time}.}
\item{data}{a data frame (must be a balanced panel).}
\item{Y}{outcome.}
\item{D}{treatment.}
\item{X}{covariates. If a covariate is time-varying, it will be averaged based on \code{X.avg.time} before balancing.}
\item{X.avg.time}{a list of time periods over which covariates are being averaged. Ignored if the treatment starts at different times.}
\item{index}{a two-element string vector specifying the unit (group)
    and time indicators. Must be of length 2.}
\item{trim.npre}{a numeric value indicating the smallest number of pre-treatment periods for a treated unit to be preserved. The default is 0.}
\item{Y.match.time}{a set of pre-treatment time periods in which the outcome variable is being balanced on.}
\item{Y.match.npre}{a numeric value indicating the number of pre-treatment outcome periods to be balanced on. If \code{Y.match.npre = 0}, no pre-treatment outcome will be part of the balancing scheme.}
\item{demean}{a logical flag indicating whether a demeaning procedure will be performed to take out the average of pre-treatment outcomes for each unit.}
\item{estimator}{a string specifying the balancing approach: "mean" for mean-balancing or "kernel" (default) for kernel-balancing via KMMD minimization.}
\item{sigma}{a numeric value specifying the bandwidth of the Gaussian kernel. If NULL (default), the median heuristic is used. Ignored if \code{estimator = "mean"}.}
\item{eta}{entropy regularization parameter for the KMMD solver. Default is 0.01.}
\item{nystrom}{logical or NULL. If TRUE, use Nystrom approximation for the kernel matrix. If NULL (default), automatically enabled when Nco > 500.}
\item{nystrom_m}{integer specifying the number of landmark points for Nystrom approximation. Default is 200.}
\item{print.baltable}{a logical flag that controls whether to print out the balance table after the algorithm is run. Ignored if treatment timing is different.}
\item{vce}{a string specifying the variance estimator. \code{vce = "none"}: no uncertainty estimates; \code{vce = "fixed.weights"}: treating balancing weights as fixed; \code{vce = "bootstrap"}: conducting non-parametric bootstrapping; \code{vce = "jackknife"}: conducting jackknife by omitting one treated unit at a time. When the treatment timing is different, \code{vce = "jackknife"} is the only available choice.}
\item{conf.lvl}{a positive number in the range of 0 and 1 specifying confidence levels for uncertainty estimates. The default value is 0.95.}
\item{nsims}{an integer specifying the number of bootstrap runs. Default is 200. Ignored if \code{vce = "none"} or \code{vce = "jackknife"}.}
\item{parallel}{a logical flag indicating whether parallel computing will be used in bootstrap/jackknife simulations.}
\item{cores}{an integer indicating the number of cores to be used in parallel computing. If not specified, the algorithm will use the maximum number of logical cores available.}
\item{seed}{an integer that sets the seed in random number generation. Ignored if \code{vce = "none"}.}
\item{\dots}{additional arguments.}
}
\details{
  \code{tjbal} provides a general reweighting approach to causal inference with time-series cross-sectional (TSCS) data. It reweights control units by minimizing the kernel Maximum Mean Discrepancy (KMMD) between treated and control groups with entropy regularization. This seeks approximate balance on a kernel-based feature expansion of the pre-treatment outcomes and covariates. The approach inherits the ability of synthetic control and latent factor models to tolerate time-varying confounders, but (1) improves feasibility and stability with reduced user discretion; (2) accommodates both short and long pre-treatment time periods with many or few treated units; and (3) balances on the high-order "trajectory" of pre-treatment outcomes rather than their period-wise average.
}
\value{
  \item{att}{average treatment effect on the treated over time.}
  \item{att.avg}{average treatment effect on the treated (averaged both across unit and over time).}
  \item{weights.co}{weights of the control units; sum to 1.}
  \item{est.att}{inference for \code{att}.}
  \item{est.att.avg}{inference for \code{att.avg}.}
  \item{bal.table}{balance table.}
  \item{data.wide}{a matrix storing data in wide form.}
  \item{id.tr}{a vector of row numbers for the treatment units.}
  \item{id.co}{a vector of row numbers for the control units.}
  \item{Y.tr}{data of the treated unit outcome.}
  \item{Y.co}{data of the control unit outcome.}
  \item{K}{the kernel matrix used for balancing.}
}
\author{
  Chad Hazlett <chazlett@ucla.edu>, UCLA

  Yiqing Xu <yiqingxu@stanford.edu>, Stanford

  Samuel G.Z. Asher (Maintainer) <sgzasher@stanford.edu>, Stanford
}
\references{
  Hazlett, Chad and Yiqing Xu, 2020. "Trajectory Balancing: A General Reweighting Approach to Causal Inference with Time-Series Cross-Sectional Data." Working Paper, UCLA and Stanford.
}
\seealso{
  \code{\link{plot.tjbal}} and \code{\link{print.tjbal}}
}
\examples{
library(fast.tjbal)
data(tjbal)
out <- tjbal(roa ~ treat + so_portion + rev2007, data = npc,
       index = c("gvkey","fyear"), estimator = "mean")
print(out)
}
