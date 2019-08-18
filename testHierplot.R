
Hierplot2<-function(hierdata, ptnum=10, param="risk difference", OR_xlim=c(0,5) ){
  # hierdata is the result from function Hier
  # ptnum is the number of AE we want to plot
  # param is the summary statistic we use to select the AE, it can be either "risk difference" or "odds ratio"
  # for this function, it will first select the top ptnum number of AE based on the summary statistic param from the selected method
  # then it will plotted the mean, 2.5% quantile, 97.5% quantile of the param of these AE
  # OR_ylim is for user to set up the y-axis limit for plotting based on "odds ratio"


  library(ggplot2)
  library(data.table)
  inputdata<-copy(hierdata)

  if (param=="risk difference"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "Diff_mean"], decreasing=TRUE), ]
    test<-head(test, ptnum)

    # change the name for plotting
    setnames(test, old=c("Diff_2.5%","Diff_97.5%" ), new=c("Diff_L","Diff_U"))

    # create two set of test
    test<-rbind(test, test)

    # create the y column for plotting
    test$yloc <- c(seq(from=ptnum+0.15, to=1.15, by=-1), seq(ptnum-0.15, to=0.85, by=-1))

    # create new column New_Diff, New_Diff_L, New_Diff_U
    # with the data from model in above and data from Raw in bottom
    test$New_Diff[1:ptnum]<-test$Diff_mean[1:ptnum]
    test$New_Diff[(ptnum+1):(2*ptnum)]<-test$Raw_Risk_Diff[(ptnum+1):(2*ptnum)]
    test$New_Diff_L[1:ptnum]<-test$Diff_L[1:ptnum]
    test$New_Diff_L[(ptnum+1):(2*ptnum)]<-NA
    test$New_Diff_U[1:ptnum]<-test$Diff_U[1:ptnum]
    test$New_Diff_U[(ptnum+1):(2*ptnum)]<-NA
    test$group<-c(rep("Model based", ptnum), rep("Raw data", ptnum))
    test$textAE<-paste0(test$AEt, "/", test$Nt, " VS ", test$AEc, "/", test$Nc)
    test$textAE[1:ptnum]<-NA

    library(ggplot2)
    # add the points at Mean
    p<-ggplot(test, aes(x=New_Diff, y=yloc, shape=group))+geom_point(size=2)
    # add line to shown confidence interval
    p<-p+geom_segment(aes(x=New_Diff_L, xend=New_Diff_U, y=yloc, yend=yloc), size=1, na.rm=TRUE)

    # add number of occurence on
    p<-p+geom_text(aes(label=textAE, x=New_Diff+0.03, y=yloc), size=3, show.legend = FALSE, na.rm = TRUE)

    # add title
    p<-p+ggtitle(paste0("Top ", ptnum, " AE of mean risk difference plotted with 95% credible interval"))
    p<-p + theme(plot.title = element_text(size=15, hjust=0.5))

    # ylable and x label
    p<-p+xlab("Risk Difference")+ylab("PT")
    p<-p + theme(axis.title.x = element_text(size=13)) + theme(axis.title.y = element_text(size=13))

    # x-axis coordinate and y-axis coordinate
    p<-p + theme(axis.text.y = element_text(color = as.factor(test$SoC[1:ptnum]), size = 13))
    p<-p + scale_y_continuous(breaks=seq(from=ptnum,to=1,by=-1),labels=test$PT[1:ptnum])
    p<-p + theme(axis.text.x=element_text(size=13))

    # legend size
    p<-p+theme(legend.text = element_text(size=15), legend.title = element_blank())

    return(p)

  }

  if (param=="odds ratio"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "OR_median"], decreasing=TRUE), ]
    test<-head(test, ptnum)

    # change the name for plotting
    setnames(test, old=c("OR_2.5%","OR_97.5%" ), new=c("OR_L","OR_U"))

    # create two set of test
    test<-rbind(test, test)

    # create the y column for plotting
    test$yloc <- c(seq(from=ptnum+0.15, to=1.15, by=-1), seq(ptnum-0.15, to=0.85, by=-1))

    # create new column New_Diff, New_Diff_L, New_Diff_U
    # with the data from model in above and data from Raw in bottom
    test$New_OR[1:ptnum]<-test$OR_median[1:ptnum]
    test$New_OR[(ptnum+1):(2*ptnum)]<-test$Raw_OR[(ptnum+1):(2*ptnum)]
    test$New_OR_L[1:ptnum]<-test$OR_L[1:ptnum]
    test$New_OR_L[(ptnum+1):(2*ptnum)]<-NA
    test$New_OR_U[1:ptnum]<-test$OR_U[1:ptnum]
    test$New_OR_U[(ptnum+1):(2*ptnum)]<-NA
    test$group<-c(rep("Model based", ptnum), rep("Raw data", ptnum))
    test$textAE<-paste0(test$AEt, "/", test$Nt, " VS ", test$AEc, "/", test$Nc)
    test$textAE[1:ptnum]<-NA
    # remove the Raw_OR that is Inf or out of the bounder of OR_xlim
    for(i in (ptnum+1):(2*ptnum)){
      if (test$New_OR[i]>OR_xlim[2]){
        test$textAE[i-ptnum]<-NA
        test$New_OR[i]<-NA
      }
    }

    TEST<-copy(test)

    library(ggplot2)
    # add the points at Mean
    p<-ggplot(test, aes(x=New_OR, y=yloc, shape=group))+geom_point(size=2, na.rm=TRUE) + coord_cartesian(xlim = OR_xlim)
    # add line to shown confidence interval
    p<-p+geom_segment(aes(x=New_OR_L, xend=New_OR_U, y=yloc, yend=yloc), size=1, na.rm=TRUE)

    # add number of occurence on
    p<-p+geom_text(aes(label=textAE, x=New_OR+((OR_xlim[2]-OR_xlim[1])/10), y=yloc), size=3, show.legend = FALSE, na.rm = TRUE)

    # add title
    p<-p+ggtitle(paste0("Top ", ptnum, " AE of median odds ratio plotted with 95% credible interval"))
    p<-p + theme(plot.title = element_text(size=15, hjust=0.5))

    # ylable and x label
    p<-p+xlab("Risk Difference")+ylab("PT")
    p<-p + theme(axis.title.x = element_text(size=13)) + theme(axis.title.y = element_text(size=13))

    # x-axis coordinate and y-axis coordinate
    p<-p + theme(axis.text.y = element_text(color = as.factor(test$SoC[1:ptnum]), size = 13))
    p<-p + scale_y_continuous(breaks=seq(from=ptnum,to=1,by=-1),labels=test$PT[1:ptnum])
    p<-p + theme(axis.text.x=element_text(size=13))

    # legend size
    p<-p+theme(legend.text = element_text(size=15), legend.title = element_blank())

    return(p)
  }
}
