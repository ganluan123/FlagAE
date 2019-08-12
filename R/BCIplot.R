#########################################################################################################
#########################################################################################################

#' @name BCIplot
#'
#' @title Plot based on Fisher exact test
#'
#' @description Both functions \code{BCIplot} plot the top AEs with
#' top ptnum (one parameter, an integer) highest difference of incidence rate
#' between treatment and control group(treatment - control). The incidence and
#' confidence interval from treatment and control group are plotted. AEs from
#' the same SoC is plotted in same color. Function \code{BCItable} return a table containing
#' information about AE selected.
#'
#' @details The incidence and confidence interval for both control and treatment
#' group are calculated by Fisher exact test. \code{BCIplot} plots out the AEs with confidence
#' interval, while \code{BCItable} output a table containig the information of
#' AEs in the plot. The table from \code{BCItable} has the following columns:\
#' AEBODSYS: SOC of the AE
#' AEDECOD: PT term of the AE
#' AEc: number of occurence of the AE in control group
#' AEt: number of occurence of the AE in treatment group
#' Nc: number of subjects in control group
#' Nt: number of subjects in treatment group
#' c_mean: incidence of the AE in control group
#' c_lowerbd: lower bound of the confidence interval of incidence the AE in control group
#' c_upperbd: upper bound of the confidence interval of incidence the AE in control group
#' t_mean: incidence of the AE in treatment group
#' t_lowerbd: lower bound of the confidence interval of incidence the AE in treatment group
#' t_upperbd: upper bound of the confidence interval of incidence the AE in treatment group
#' mean_diff: difference of incidence of the AE in treatment group and control group
#'
#'
#' @param aedata output from function \code{\link{preprocess}}
#' @param ptnum positive integer, number of AEs to be plotted, default is 10
#' @param conf.level number between 0 and 1, confidence level for constructing
#'   confidence interval, default is 0.95
#'
#'
#' @examples
#' \dontrun{
#' data(ADAE)
#' data(ADSL)
#' AEdata<-preprocess(adsl=ADSL, adae=ADAE)
#' BCIplot(aedata=AEdata)
#' BCIplot(aedata=AEdata, ptnum=15, conf.level=0.9)
#' BCItable(aedata=AEdata)
#' BCItable(aedata=AEdata, ptnum=15, conf.level=0.9)
#' }
#' @seealso \code{\link{preprocess}}
#'
#' @export
BCIplot <- function(aedata, ptnum=10, conf.level=0.95) {

  # This function calculate the incidence rate for each adverse event (AE) in
  # treatment and in control group and also the confidence interval for this
  # incidence getting from fisher exact test then the AEs with the highest ptnum
  # (with default be 10) incidence rate difference between treatment and control
  # group were plotted with incidence and conf. interval shown pt from the same
  # SoC group were plotted in the same color.

  # aedata: the output from function preprocess
  # ptnum: the number of top AEs we want to plot, with default be 10
  # conf.level: confident level used for fisher exact test with default be 0.95

  library(binom)
  # the incidence and low, up bound of CI for control group
  aedata$mc <- aedata$AEc/aedata$Nc
  aedata$llc <- binom.confint(aedata$AEc,aedata$Nc, conf.level, method="exact")$lower
  aedata$ulc <- binom.confint(aedata$AEc,aedata$Nc, conf.level, method="exact")$upper

  # the incidence and low, up bound of CI for treatment group
  aedata$mt <- aedata$AEt/aedata$Nt
  aedata$llt <- binom.confint(aedata$AEt,aedata$Nt, conf.level, method="exact")$lower
  aedata$ult <- binom.confint(aedata$AEt,aedata$Nt, conf.level, method="exact")$upper

  # get the difference of incidence between control and treatment group
  aedata$d <- abs(aedata$mt - aedata$mc)
  # order the dataset by the difference in descending order
  aedata1 <- aedata[order(-aedata$d),]

  # change the datatype to character
  pdat <- data.frame(lapply(aedata1, as.character), stringsAsFactors=FALSE)

  # start the plotting
  plot(NA, xlim=c(-1.1,1.1), ylim=c(-5,5*ptnum), main=paste0("Confidence interval of top ", ptnum," AE in incidene difference by treatment group"),
       xlab="AE rate", ylab="PT", yaxt='n')

  # plot PT from same SOC in same color
  Soc1<-unique(pdat[1:ptnum, "b"])

  for (iae in 1:ptnum){
    # find the color for the PT
    for (j in 1:length(Soc1)){
      if (pdat[iae, "b"]==Soc1[j]) Col<-j
    }

    # determine the position of these plotted AE on y-axis
    # with AE of the highest incidence difference between treatment and control on the top
    yc<-5*ptnum-5*iae
    yt<-yc-1

    lines(c(pdat[iae,]$llc,pdat[iae,]$ulc),c(yc,yc))
    points(pdat[iae,]$mc,yc,pch=4,cex = 0.6)
    lines(c(pdat[iae,]$llt,pdat[iae,]$ult),c(yt,yt))
    points(pdat[iae,]$mt,yt,cex = 0.6,pch=16)
    # text(0.02, (yc+yt)/2, pdat[iae,]$b,cex = 0.9,adj=1)
    text(-0.05, (yc+yt)/2, pdat[iae,]$AEDECOD,cex = 0.9,adj=1, col=Col)
  }

  legend(0.6, 3.5*ptnum, legend=c("control", "treatment"), pch=c(4, 16), lty=c(1,1))

}


#########################################################################################################
#########################################################################################################
#' @rdname BCIplot
#' @export

BCItable<-function(aedata, ptnum=10, conf.level=0.95){
  # this fucntion take the same parameter as function BCIplot and
  # does the same thing as BCIplot, instead of plotting out the
  # AEs this function return a table containg the detailed information
  # about this top AEs

  library(binom)
  # the incidence and low, up bound of CI for control group
  aedata$c_mean <- aedata$AEc/aedata$Nc
  aedata$c_lowerbd <- binom.confint(aedata$AEc,aedata$Nc, conf.level, method="exact")$lower
  aedata$c_upperbd <- binom.confint(aedata$AEc,aedata$Nc, conf.level, method="exact")$upper

  # the incidence and low, up bound of CI for treatment group
  aedata$t_mean <- aedata$AEt/aedata$Nt
  aedata$t_lowerbd <- binom.confint(aedata$AEt,aedata$Nt, conf.level, method="exact")$lower
  aedata$t_upperbd <- binom.confint(aedata$AEt,aedata$Nt, conf.level, method="exact")$upper

  # get the difference of incidence between control and treatment group
  aedata$mean_diff <- abs(aedata$t_mean - aedata$c_mean)
  # order the dataset by the difference in descending order
  aedata1 <- aedata[order(-aedata$mean_diff),]

  # change the datatype to character
  pdat <- data.frame(lapply(aedata1, as.character), stringsAsFactors=FALSE)

  # drop the column b, i, j
  drops<-c("b", "i", "j")
  pdat<-pdat[, !names(pdat) %in% drops]

  # return the AE with top ptnum mean difference
  head(pdat, ptnum)
}



