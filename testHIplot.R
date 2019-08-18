HIplot2<-function(hierdata,isingdata, aedata, ptnum=10, param="risk difference", OR_xlim=c(0,5) ){
  # hierdata is the result from function Hier
  # isingdata is the result from function Ising
  # aedata is the result from function preprocess
  # ptnum is the number of AE we want to plot
  # param is the summary statistic we use to select the AE, it can be either "risk difference" or "odds ratio"
  # for this function, it will first select the top ptnum number of AE based on the summary statistic param from the both models(methods)
  # then it will plotted the mean, 2.5% quantile, 97.5% quantile of the param of these AE from both two models(methods)
  # it also has one more parameter aedata
  # this is for fucntion BCItable to get the top AEs with risk difference based on fisher exact test
  # OR_ylim is for user to set up the y-axis limit for plotting based on "odds ratio"

  library(ggplot2)
  library(data.table)

  hierdata.new<-copy(hierdata)
  hierdata.new$gammaeq1<-NA
  inputdata<- rbind(hierdata.new,isingdata)

  # get the PTs of AE with top ptnum difference in risk difference based on fisher exact test
  BCI<-BCItable(aedata, ptnum)$AEDECOD

  if (param=="risk difference"){
    # first to get the top 10 AEs from each method
    test1_Hier<-subset(inputdata, Method=="Bayesian Hierarchical Model")
    test2_Hier<-test1_Hier[order(test1_Hier[, "Diff_mean"], decreasing=TRUE), ]
    test3_Hier<-head(test2_Hier, ptnum)

    test1_Is<-subset(inputdata, Method=="Ising prior")
    test2_Is<-test1_Is[order(test1_Is[, "Diff_mean"], decreasing=TRUE), ]
    test3_Is<-head(test2_Is, ptnum)

    # get the PT selected by each method
    PT_Hier<-test3_Hier$PT
    PT_Is<-test3_Is$PT

    # get the union of PTs, and PTs selected unique by each method
    PT_Intsect<-intersect(PT_Hier, PT_Is)
    PT_Hier_Uni<-setdiff(PT_Hier, PT_Intsect)
    PT_Is_Uni<-setdiff(PT_Is, PT_Intsect)

    # get the number of PTs in each set
    n_Intsect<-length(PT_Intsect)
    n_Hier_Uni<-length(PT_Hier_Uni)
    n_Is_Uni<-length(PT_Is_Uni)

    # get the indicator for PT in different set
    # 0 for PT in union, 1 for PT unique in Hierarchical model
    # 2 for PT in Ising model
    test_Intsect<-inputdata[which(inputdata$PT %in% PT_Intsect),]
    test_Intsect$Ind<-"A"
    test_Intsect$Set<-"Intersect of both models"
    test_Hier_Uni<-inputdata[which(inputdata$PT %in% PT_Hier_Uni),]
    test_Hier_Uni$Ind<-"B"
    test_Hier_Uni$Set<-"Unique in Hierarchical Model"
    test_Is_Uni<-inputdata[which(inputdata$PT %in% PT_Is_Uni),]
    test_Is_Uni$Ind<-"C"
    test_Is_Uni$Set<-"Unique in Ising prior"

    # combine the dataset together
    test<-rbind(test_Intsect, test_Hier_Uni)
    test<-rbind(test, test_Is_Uni)

    # add one more column "category" for the indication of whether the AE is inside BCI, AE slected by fisher exact test
    test$category<-0
    if(length(intersect(BCI, test$PT)))  test[which(test$PT %in% BCI), ]$category<-1
    # the if statement is to make sure there is at least one PT that in both BCI and test

    # sort the test by Method, Ind, and -Diff_mean
    test<-test[order(test$Method, test$Ind, -test$Diff_mean),]
    # add half of test to the end
    temptest<-copy(test[1:(dim(test)[1]/2),])
    temptest$Method<-"Raw Data"
    temptest$textAE<-paste0(temptest$AEt, "/", temptest$Nt, " VS ", temptest$AEc, "/", temptest$Nc)
    # temptest$Diff_mean<-NA
    temptest$`Diff_2.5%`<-NA
    temptest$`Diff_97.5%`<-NA
    # delete information of Raw_Risk_Diff
    # test$Raw_Risk_Diff<-NA
    test$textAE<-NA

    # combine two sets together
    test<-rbind(test, temptest)
    # get the number of sets of plot
    N<-dim(test)[1]/3
    # create the new variable
    test$New_Diff<-c(test$Diff_mean[1:(2*N)], test$Raw_Risk_Diff[(2*N+1):(3*N)])
    # change the name for plotting
    setnames(test, old=c("Diff_2.5%","Diff_97.5%" ), new=c("Diff_L","Diff_U"))



    # create yloc for the location of the plot
    test$yloc<-c(seq(from=1.5*N+0.3, to=1.5+0.3, by=-1.5), seq(from=1.5*N+0.15, to=1.5+0.15, by=-1.5),
                 seq(from=1.5*N, to=1.5, by=-1.5))

    library(ggplot2)
    # add the points at Mean
    p<-ggplot(test, aes(x=New_Diff, y=yloc, shape=Method, color=Set))+geom_point(size=2)
    # add line to shown confidence interval
    p<-p+geom_segment(aes(x=Diff_L, xend=Diff_U, y=yloc, yend=yloc, color=Set), size=1, na.rm=TRUE)

    # add number of occurence on
    p<-p+geom_text(aes(label=textAE, x=New_Diff+0.03, y=yloc), size=3, show.legend = FALSE, na.rm = TRUE)

    # add title
    p<-p+ggtitle(paste0("Top ", ptnum, " AE of mean risk difference "))
    p<-p + theme(plot.title = element_text(size=15, hjust=0.5))

    # ylable and x label
    p<-p+xlab("Risk Difference")+ylab("PT")
    p<-p + theme(axis.title.x = element_text(size=13)) + theme(axis.title.y = element_text(size=13))

    # x-axis coordinate and y-axis coordinate
    p<-p + theme(axis.text.y = element_text(color = as.factor(test$category), size = 13))
    # p<-p + theme(axis.text.y = element_text(color = as.factor(test$PT[1:N]), size = 13))
    p<-p + scale_y_continuous(breaks=seq(from=N*1.5,to=1.5,by=-1.5),labels=test$PT[1:N])
    p<-p + theme(axis.text.x=element_text(size=13))

    # legend size
    p<-p+theme(legend.text = element_text(size=15), legend.title = element_blank())

    return(p)
  }

  if (param=="odds ratio"){
    # first to get the top 10 AEs from each method
    test1_Hier<-subset(inputdata, Method=="Bayesian Hierarchical Model")
    test2_Hier<-test1_Hier[order(test1_Hier[, "OR_median"], decreasing=TRUE), ]
    test3_Hier<-head(test2_Hier, ptnum)

    test1_Is<-subset(inputdata, Method=="Ising prior")
    test2_Is<-test1_Is[order(test1_Is[, "OR_median"], decreasing=TRUE), ]
    test3_Is<-head(test2_Is, ptnum)

    # get the PT selected by each method
    PT_Hier<-test3_Hier$PT
    PT_Is<-test3_Is$PT

    # get the union of PTs, and PTs selected unique by each method
    PT_Intsect<-intersect(PT_Hier, PT_Is)
    PT_Hier_Uni<-setdiff(PT_Hier, PT_Intsect)
    PT_Is_Uni<-setdiff(PT_Is, PT_Intsect)

    # get the number of PTs in each set
    n_Intsect<-length(PT_Intsect)
    n_Hier_Uni<-length(PT_Hier_Uni)
    n_Is_Uni<-length(PT_Is_Uni)

    # get the indicator for PT in different set
    # 0 for PT in union, 1 for PT unique in Hierarchical model
    # 2 for PT in Ising model
    test_Intsect<-inputdata[which(inputdata$PT %in% PT_Intsect),]
    test_Intsect$Ind<-"A"
    test_Intsect$Set<-"Intersect of both models"
    test_Hier_Uni<-inputdata[which(inputdata$PT %in% PT_Hier_Uni),]
    test_Hier_Uni$Ind<-"B"
    test_Hier_Uni$Set<-"Unique in Hierarchical Model"
    test_Is_Uni<-inputdata[which(inputdata$PT %in% PT_Is_Uni),]
    test_Is_Uni$Ind<-"C"
    test_Is_Uni$Set<-"Unique in Ising prior"

    # combine the dataset together
    test<-rbind(test_Intsect, test_Hier_Uni)
    test<-rbind(test, test_Is_Uni)

    # add one more column "category" for the indication of whether the AE is inside BCI, AE slected by fisher exact test
    test$category<-0
    if(length(intersect(BCI, test$PT)))  test[which(test$PT %in% BCI), ]$category<-1
    # the if statement is to make sure there is at least one PT that in both BCI and test

    # sort the test by Method, Ind, and -Diff_mean
    test<-test[order(test$Method, test$Ind, -test$OR_median),]

    # add half of test to the end
    temptest<-copy(test[1:(dim(test)[1]/2),])
    temptest$Method<-"Raw Data"
    temptest$textAE<-paste0(temptest$AEt, "/", temptest$Nt, " VS ", temptest$AEc, "/", temptest$Nc)
    # temptest$Diff_mean<-NA
    temptest$`OR_2.5%`<-NA
    temptest$`OR_97.5%`<-NA
    # delete information of Raw_Risk_Diff
    # test$Raw_Risk_Diff<-NA
    test$textAE<-NA

    # combine two sets together
    test<-rbind(test, temptest)
    # get the number of sets of plot
    N<-dim(test)[1]/3
    # create the new variable
    test$New_OR<-c(test$OR_median[1:(2*N)], test$Raw_OR[(2*N+1):(3*N)])
    # change the name for plotting
    setnames(test, old=c("OR_2.5%","OR_97.5%" ), new=c("OR_L","OR_U"))

    # create yloc for the location of the plot
    test$yloc<-c(seq(from=1.5*N+0.3, to=1.5+0.3, by=-1.5), seq(from=1.5*N+0.15, to=1.5+0.15, by=-1.5),
                 seq(from=1.5*N, to=1.5, by=-1.5))

    # remove the Raw_OR that is Inf or out of the bounder of OR_xlim
    for(i in (2*N+1):(3*N)){
      if (test$New_OR[i]>OR_xlim[2]){
        test$textAE[i]<-NA
        test$New_OR[i]<-NA
      }
    }

    library(ggplot2)
    # add the points at Mean
    p<-ggplot(test, aes(x=New_OR, y=yloc, shape=Method, color=Set))+geom_point(size=2, na.rm=TRUE)+ coord_cartesian(xlim = OR_xlim)
    # add line to shown confidence interval
    p<-p+geom_segment(aes(x=OR_L, xend=OR_U, y=yloc, yend=yloc, color=Set), size=1, na.rm=TRUE)

    # add number of occurence on
    p<-p+geom_text(aes(label=textAE, x=New_OR+((OR_xlim[2]-OR_xlim[1])/10), y=yloc), size=3, show.legend = FALSE, na.rm = TRUE)

    # add title
    p<-p+ggtitle(paste0("Top ", ptnum, " AE of median odds ratio "))
    p<-p + theme(plot.title = element_text(size=15, hjust=0.5))

    # ylable and x label
    p<-p+xlab("Risk Difference")+ylab("PT")
    p<-p + theme(axis.title.x = element_text(size=13)) + theme(axis.title.y = element_text(size=13))

    # x-axis coordinate and y-axis coordinate
    p<-p + theme(axis.text.y = element_text(color = as.factor(test$category), size = 13))
    # p<-p + theme(axis.text.y = element_text(color = as.factor(test$PT[1:N]), size = 13))
    p<-p + scale_y_continuous(breaks=seq(from=N*1.5,to=1.5,by=-1.5),labels=test$PT[1:N])
    p<-p + theme(axis.text.x=element_text(size=13))

    # legend size
    p<-p+theme(legend.text = element_text(size=15), legend.title = element_blank())

    return(p)
  }
}
