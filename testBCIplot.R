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
  aedata1 <- copy(aedata[order(-aedata$d),])

  # change the datatype to character
  pdat <- data.frame(lapply(aedata1, as.character), stringsAsFactors=FALSE)
  
  # get the top ptnum AE by risk difference 
  pdat <- head(pdat, ptnum)
  
  # stack two pdat together
  pdat<-rbind(pdat, pdat)
  
  # give y values to the dataset
  pdat$yloc <- c(seq(from=ptnum+0.15, to=1.15, by=-1), seq(ptnum-0.15, to=0.85, by=-1))
  
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



