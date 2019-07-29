#' @name Data
#'
#' @title Example dataset
#'
#' @description
#' Randomly generated datasets used for running examples in
#' this package.
#'
#' @usage data(ADAE)
#'
#' @details \code{ADAE} is the "adverse event analysis dataset", it contains four columns:\cr
#' \code{USUBJID}: unique subject id \cr
#' \code{AEBODSYS}: SoC of the AE \cr
#' \code{AEDECOD}: PT of the AE \cr
#' \code{TRTEMFL}: indicitator of whether adverse events occurs within 30 days after taking the first dose \cr
#'
#' \code{ADSL} is the "subject level analysis dataset", it contains three columns: \cr
#' \code{USUBJID}: unique subject id \cr
#' \code{TREATMENT}: treatment applied on each subject, subjects are divided into control or treatment group based on this
#' \code{SAFFL}: indicitator that patients at least took one dose of the treatment
#'
#' @docType data
#'
#' @keywords datasets
#'
"ADAE"


#' @rdname Data
#'
#' @usage data(ADSL)
"ADSL"
