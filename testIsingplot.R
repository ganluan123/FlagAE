Isingplot<-function(isingdata, ptnum=10, param="risk difference" ){
  # isingdata is the result from function Ising
  # ptnum is the number of AE we want to plot
  # param is the summary statistic we use to select the AE, it can be either "risk difference" or "odds ratio"
  # for this function, it will first select the top ptnum number of AE based on the summary statistic param from the selected method
  # then it will plotted the mean, 2.5% quantile, 97.5% quantile of the param of these AE


  library(ggplot2)
  library(data.table)
  inputdata<-isingdata

  if (param=="risk difference"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "Diff_mean"], decreasing=TRUE), ]
    test<-head(test, ptnum)

    # create the x column for plotting
    len<-dim(test)[1]
    PT_lable<-test$PT[1:len]
    test$x<-seq(1, len, by=1)

    # plot PT from same SOC in same color
    test$Color<-NA
    Soc1<-unique(test$b)

    for (iae in 1:ptnum){
      # find the color for the PT
      for (j in 1:length(Soc1)){
        if (test[iae, "b"]==Soc1[j]) test[iae, "Color"]<-j
      }
    }

    # change the name for plotting
    setnames(test, old=c("Diff_2.5%","Diff_97.5%" ), new=c("Diff_L","Diff_U"))

    p <- ggplot(test, aes(x=x, y=Diff_mean)) + geom_pointrange(aes(ymin=Diff_L, ymax=Diff_U))
    p1 <-p + labs(x = "Prefered Term",y="Risk difference",title = paste0("Top ", ptnum, " AE of mean risk difference plotted with 95% interval"))


    p2 <- p1 + theme(axis.text.x = element_text(color = test$Color, size = 8, angle = 60, hjust = 1), axis.text.y = element_text(color = "black", size = 8))
    p3 <- p2 + scale_x_continuous(breaks=seq(1,len,by=1),labels=PT_lable)
    return(p3)
  }

  if (param=="odds ratio"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "OR_mean"], decreasing=TRUE), ]
    test<-head(test, ptnum)

    # create the x column for plotting
    len<-dim(test)[1]
    PT_lable<-test$PT[1:len]
    test$x<-seq(1, len, by=1)

    # plot PT from same SOC in same color
    test$Color<-NA
    Soc1<-unique(test$b)

    for (iae in 1:ptnum){
      # find the color for the PT
      for (j in 1:length(Soc1)){
        if (test[iae, "b"]==Soc1[j]) test[iae, "Color"]<-j
      }
    }

    # change the name for plotting
    setnames(test, old=c("OR_2.5%", "OR_97.5%" ), new=c("OR_L", "OR_U"))

    p <- ggplot(test, aes(x=x, y=OR_mean)) + geom_pointrange(aes(ymin=OR_L, ymax=OR_U)) + coord_cartesian(ylim = c(0, 5))
    p1 <-p + labs(x = "Prefered Term",y="Odds ratio",title = paste0("Top ", ptnum, " AE of mean odds ratio plotted with 95% interval"))


    p2 <- p1 + theme(axis.text.x = element_text(color = test$Color, size = 8, angle = 60, hjust = 1), axis.text.y = element_text(color = "black", size = 8))
    p3 <- p2 + scale_x_continuous(breaks=seq(1,len,by=1),labels=PT_lable)
    return(p3)

  }
}

Isingtable<-function(isingdata, ptnum=10, param="risk difference" ){
  # this function takes the same parameter as function Hierplot
  # it returns the detailed information of AEs plotted in Hierplot

  library(data.table)
  inputdata<-isingdata

  if (param=="risk difference"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "Diff_mean"], decreasing=TRUE), ]
    test<-head(test, ptnum)
  }

  if (param=="odds ratio"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "OR_mean"], decreasing=TRUE), ]
    test<-head(test, ptnum)
  }
  return(test)
}
