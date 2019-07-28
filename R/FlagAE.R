#' FlagAE: A package for flagging adverse events by Bayesian methods.
#'
#' @details This packages can flag adverse events by three stage Bayesian Hierarchical model (see \code{\link{Hier}}) and by
#' Bayesian model with Ising latent variable (see \code{\link{Ising}}). It also provide the option to flag adverase events by
#' frequentist method (see \code{\link{gci2}}). Adverse events selected by both Bayesian methods can be plotted together
#' for comparison (see \code{\link{HIPLOT}}). Last this package provides the options to compare these two bayesian models
#' or compare differnt prior/initial values/parameter values in one model by calculating the loss function with
#' cross validation (see \code{\link{CVhier}}, \code{\link{CVising}}).\cr
#' Two datasets are also provided in this packages for demo purpose. They are \code{\link{ADAE}} (adverse event analysis dataset),
#' and \code{\link{ADSL}} (subject level analysis dataset). They can be used for running the example code for each function to give
#' users a better understanding of these functions.\cr
#'
#' @section Functions in this package:
#'
#' Raw data process: \code{\link{preprocess}}, \code{\link{preprocess2}}\cr
#' Plot based on Fisher exact test: \code{\link{gci}}, \code{\link{gci2}} \cr
#' Three stages Bayesian Hierarchical Model: \code{\link{Hier_history}}, \code{\link{sum_Hier}}, \code{\link{Hier}}, \code{\link{Hiergetpi}} \cr
#' Bayesian Model with Ising Latent Variables: \code{\link{Ising_history}}, \code{\link{sum_Ising}}, \code{\link{Ising}}, \code{\link{Isinggetpi}} \cr
#' Plot for two Bayesian models: \code{\link{gci3}}, \code{\link{HIPLOT}} \cr
#' Loss Calculation by Cross Validation: \code{\link{Lossfun}}, \code{\link{kfdpar}}, \code{\link{CVhier}}, \code{\link{cvising}}\cr
#'
#' @section Dataset in this package:
#'
#' \code{\link{ADAE}}, \code{\link{ADSL}}
#'
#' @section preinstalled packages required:
#'
#' dplyr, tidyr, binom, ggplot2, data.table, foreach, doParallel,
#' mcmcplots, rjags, R2jags
#'
#'
#'
#' @docType package
#'
#' @name FlagAE
NULL
