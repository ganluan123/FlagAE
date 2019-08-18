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
  
  # get the top ptnum AE by risk difference 
  pdat <- head(pdat, ptnum)
  
  # stack two pdat together
  pdat<-rbind(pdat, pdat)
  
  # give y values to the dataset
  pdat$yloc <- c(seq(from=ptnum+0.15, to=1.15, by=-1), seq(ptnum-0.15, to=0.85, by=-1))
  # give y value for text display
  pdat$ytext<-c(seq(from=ptnum, to=1, by=-1), seq(from=ptnum, to=1, by=-1))
  # delete the AEDECOD of the second half of padt for text display purpose
  pdat$AEDECOD[(ptnum+1):(2*ptnum)]<-NA
  
  # add the text for showing number of number of occurane and total number of patients
  pdat$textAE[1:ptnum]<-paste0(pdat$AEc[1:ptnum], "/", pdat$Nc[1:ptnum])
  pdat$textAE[(ptnum+1):(2*ptnum)]<-paste0(pdat$AEt[(ptnum+1):(2*ptnum)], "/", pdat$Nt[(ptnum+1):(2*ptnum)])
  
  # create new column Mean, meanlb, meanub for both control and treatment 
  # with the first ptnum rows for control and rows from ptnum+1 to 2*ptnum for treatment group
  pdat$Mean<-as.numeric(c(pdat[1:ptnum, "mc"], pdat[(ptnum+1):(2*ptnum), "mt"]))
  pdat$lb<-as.numeric(c(pdat[1:ptnum, "llc"], pdat[(ptnum+1):(2*ptnum), "llt"]))
  pdat$ub<-as.numeric(c(pdat[1:ptnum, "ulc"], pdat[(ptnum+1):(2*ptnum), "ult"]))
  pdat$group<-c(rep("control", ptnum), rep("treatment", ptnum))
  
  library(ggplot2)
  # add the points at Mean
  p<-ggplot(pdat, aes(x=Mean, y=yloc))+geom_point(size=2)
  # add line to shown confidence interval
  p<-p+geom_segment(aes(x=lb, xend=ub, y=yloc, yend=yloc, linetype=group), size=1)
  
  # add number of occurence on 
  p<-p+geom_text(aes(label=textAE, x=(ub+Mean)/2, y=yloc+0.15), size=3, show.legend = FALSE)
  
  # add title
  p<-p+ggtitle(paste0("Confidence interval of top ", ptnum, " AEs in difference by treatment group"))
  p<-p + theme(plot.title = element_text(size=15, hjust=0.5))
  
  # ylable and x label
  p<-p+xlab("AE rate")+ylab("PT")
  p<-p + theme(axis.title.x = element_text(size=13)) + theme(axis.title.y = element_text(size=13))
  
  # x-axis coordinate and y-axis coordinate  
  p<-p + theme(axis.text.y = element_text(color = as.factor(pdat$b[1:ptnum]), size = 13))
  p<-p + scale_y_continuous(breaks=seq(from=ptnum,to=1,by=-1),labels=pdat$AEDECOD[1:ptnum])
  p<-p + theme(axis.text.x=element_text(size=13))
  
  # legend size
  p<-p+theme(legend.text = element_text(size=15), legend.title = element_blank())
  
  return(p)
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



