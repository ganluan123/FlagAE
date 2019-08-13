#' FlagAE: A package for flagging adverse events by Bayesian methods.
#'
#' @details This packages can flag adverse events by three stage Bayesian Hierarchical model (see \code{\link{Hier}}) and by
#' Bayesian model with Ising latent variable (see \code{\link{Ising}}). It also provide the option to flag adverase events by
#' frequentist method (see \code{\link{BCIplot}}). Adverse events selected by both Bayesian methods can be plotted together
#' for comparison (see \code{\link{HIplot}}). Last this package provides the options to compare these two bayesian models
#' or compare differnt hyperparamters for prior distribution in one model by calculating the loss function with
#' cross validation (see \code{\link{CVhier}}, \code{\link{CVising}}).\cr
#' Two datasets are also provided in this packages for demo purpose. They are \code{\link{ADAE}} (adverse event analysis dataset),
#' and \code{\link{ADSL}} (subject level analysis dataset). They can be used for running the example code for each function to give
#' users a better understanding of these functions.\cr
#'
#' @section Functions in this package:
#'
#' Raw data process: \code{\link{preprocess}}\cr
#' Plot based on Binomial confidence interval: \code{\link{BCIplot}}, \code{\link{BCItable}} \cr
#' Three stages Bayesian Hierarchical Model: \code{\link{Hier_history}}, \code{\link{sum_Hier}}, \code{\link{Hier}}, \code{\link{Hiergetpi}}
#' , \code{\link{Hierplot}}, \code{\link{Hiertable}}\cr
#' Bayesian Model with Ising Latent Variables: \code{\link{Ising_history}}, \code{\link{sum_Ising}}, \code{\link{Ising}}, \code{\link{Isinggetpi}}
#' , \code{\link{Isingplot}}, \code{\link{Isingtable}} \cr
#' Plot for two Bayesian models: \code{\link{HIplot}}, \code{\link{HItable}} \cr
#' Loss Calculation by Cross Validation: \code{\link{Lossfun}}, \code{\link{kfdpar}}, \code{\link{CVhier}}, \code{\link{CVising}}\cr
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
