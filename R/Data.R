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
#'
#'
#' \code{ADSL} is the "subject level analysis dataset", it contains three columns: \cr
#' \code{USUBJID}: unique subject id \cr
#' \code{TREATMENT}: treatment applied on each subject, subjects are divided into control or treatment group based on this
#' \code{TRTCTR}: indicitator for treatment group (TRTCTR=1, TREATMENT = "xyz" or 'xyz2') or
#' control group (TRTCTR=0, TREATMENT='control' or 'control2')
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
