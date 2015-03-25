#' Bayesian methods for inference and model selection using logistic regression models
#'
#' Package uses MCMC methods for fitting logistic regression models in
#' a Bayesian framework. Uses reversible-jump MCMC to perform variable selection
#' and Bayesian model averaging.
#'
#' Package provides a function \code{bayesLog} for running logistic regression models in
#' a Bayesian framework. It also includes generic functions for subsetting, summarising
#' and plotting output.
#'
#' @docType package
#' @name BayesLog-package
#' @author TJ McKinley <t.mckinley@@exeter.ac.uk>
#' @importFrom Rcpp evalCpp
#' @useDynLib BayesLog
NULL

