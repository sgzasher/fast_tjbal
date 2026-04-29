\name{plot.tjbal}
\alias{plot.tjbal}
\title{Plotting}
\description{Visualizes results from trajectory balancing.}
\usage{\method{plot}{tjbal}(x, type = "gap", subgroup = NULL, xlim = NULL,
    ylim = NULL, xlab = NULL, ylab = NULL, count = TRUE,
    legendOff = FALSE, legend.pos = "bottom", legend.ncol = NULL,
    legend.labs = NULL, main = NULL, raw = "none", stat = "mean",
    trim = TRUE, trim.wtot = 0.9, theme.bw = TRUE, cex.main = NULL,
    cex.axis = NULL, cex.lab = NULL, cex.legend = NULL,
    cex.text = NULL, log.weights = FALSE, wmin = NULL,
    axis.adjust = FALSE, \dots)
}
\arguments{
  \item{x}{a \code{\link{tjbal}} object.}
  \item{type}{a string that specifies the type of the plot. Must be one
    of the following: "gap" (plotting the average treatment effect on
    the treated); "counterfactual" or "ct" (plotting predicted average Y(0)'s); "balance" (plotting balance of covariates and pre-treatment outcomes before and after reweighting); "weights" (plotting the distribution of weights of the control units).}
  \item{subgroup}{a numeric value indicating the subgroup for which the results will be plotted.}
  \item{xlim}{a two-element numeric vector specifying the range of x-axis.}
  \item{ylim}{a two-element numeric vector specifying the range of y-axis.}
  \item{xlab}{a string indicating the label of the x-axis.}
  \item{ylab}{a string indicating the label of the y-axis.}
  \item{count}{a logical flag indicating whether a bar plot of the number of treated units will be shown at the bottom of the plot. The default is TRUE.}
  \item{legendOff}{a logical flag controlling whether to show the legends.}
  \item{legend.pos}{a character controlling the position of the legends.}
  \item{legend.ncol}{a numeric value controlling the number of columns of the legends.}
  \item{legend.labs}{a vector of characters specifying the labels of the legends.}
  \item{main}{a string that controls the title of the plot.}
  \item{raw}{a string indicating whether or how raw data will be shown in the "counterfactual" plot.}
  \item{stat}{a string specifying the statistics in the balance plot, either "mean" or "sd".}
  \item{trim}{a logical flag controlling whether to trim the control group data in a counterfactual plot.}
  \item{trim.wtot}{a numeric value controlling the total sum of weights of the most heavily weighted controls to remain.}
  \item{theme.bw}{a logical flag controlling whether a black-white theme will be used.}
  \item{cex.main}{a numeric value controlling the font size of the title.}
  \item{cex.axis}{a numeric value controlling the font size of the axis labels.}
  \item{cex.lab}{a numeric value controlling the font size of the axis titles.}
  \item{cex.legend}{a numeric value controlling the font size of the legends.}
  \item{cex.text}{a numeric value controlling the font size of the text in the figure.}
  \item{log.weights}{a logical flag controlling whether to take a log transformation of the weights.}
  \item{wmin}{a numeric number specifying the minimum weight (in common logarithmic).}
  \item{axis.adjust}{a logic flag indicating whether to adjust labels on x-axis.}
  \item{\dots}{other arguments.}
}
\details{
  \code{plot.tjbal} visualizes the results of trajectory balancing.
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
  \code{\link{tjbal}} and \code{\link{print.tjbal}}
}
