#' Bayesian methods for inference using logistic regression models
#'
#' Package uses MCMC methods for fitting logistic regression models in
#' a Bayesian framework. 
#'
#' Package provides a function \code{bayesLog} for fitting binary logistic regression models in
#' a Bayesian framework. It also includes generic functions for subsetting, summarising
#' and plotting output.
#'
#' @docType package
#' @name BayesLog-package
#' @author TJ McKinley <t.mckinley@exeter.ac.uk>
#' @imports lme4
#' @imports RcppArmadillo
#' @imports coda
#' @importFrom Rcpp evalCpp
#' @useDynLib BayesLog
NULL

