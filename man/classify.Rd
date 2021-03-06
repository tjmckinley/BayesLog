% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify.R
\name{classify}
\alias{classify}
\alias{classify.bayesLog.pred}
\alias{summary.bayesLog.class}
\alias{print.bayesLog.class}
\alias{plot.bayesLog.class}
\title{Generates posterior distributions for sensitivity
specificity, PPV and NPV}
\usage{
classify(pred, obs, thresh, ...)

\method{classify}{bayesLog.pred}(pred, obs, thresh, ...)

\method{summary}{bayesLog.class}(object, ...)

\method{print}{bayesLog.class}(x, ...)

\method{plot}{bayesLog.class}(x, type = c("ind", "comp"),
  comp_type = c("all", "ROC", "PV"), ...)
}
\arguments{
\item{pred}{A \code{bayesLog.pred} object.}

\item{obs}{Vector of binary observations (of length \code{nobs})}

\item{thresh}{Vector of thresholds between 0 and 1 at which to generate
classification estimates.}

\item{\dots}{Not currently used.}

\item{object}{An object of class \code{bayesLog.class}}

\item{x}{An object of class \code{bayesLog.class}}

\item{type}{Character vector defining whether to plot
individual or comparative plots.}

\item{comp_type}{Character containing subplot to extract when doing comparative plotting.
If \code{"all"} then plots both ROC and Predictive Values plots, else one
can choose to produce either.}
}
\value{
An object of class \code{bayesLog.class}, which is a list
including a subset of elements: \code{sens}, \code{spec}, \code{ppv}, \code{npv}
which are \code{npred} x \code{nthresh} matrices, and and element \code{thresh}
which is a vector of thresholds.
}
\description{
Generates posterior distributions for sensitivity
specificity, PPV and NPV from a set of posterior predictions.
}
\section{Methods (by class)}{
\itemize{
\item \code{bayesLog.pred}: Classify method for \code{bayesLog.pred} objects

\item \code{bayesLog.class}: Summary method for \code{bayesLog.class} objects

\item \code{bayesLog.class}: Print method for \code{bayesLog.class} objects

\item \code{bayesLog.class}: Plot method for \code{bayesLog.class} objects
}}

