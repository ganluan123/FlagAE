#' @name HierIsingplot
#'
#' @title Plot for two Bayesian models
#'
#' @description These two functions are for plotting the top AEs selected by Bayesian hierarchical model
#' and Bayesian model with Ising prior.
#'
#' @details \code{gci3} calculates the incidence rate for each adverse event (AE) in treatment
#' and in control group and also the confidence interval for this incidence getting from fisher exact test
#' then the PT of AEs with the highest \code{ptnum} incidence rate difference between treatment and control
#' group are returned. \cr
#' \code{HIPLOT} first selects the top \code{ptnum} (an integer) AE based on the selected statistic (either "odds ratio" or "risk difference").
#' Then it plots the mean, 2.5% quantile, 97.5% quantile of the selected statistic of these AE from both two models(methods). It seperates
#' the AEs slected by both Bayesian methods from AEs selected by only one method. Also it indicates whether the AE selected by these two Bayesian models were also selected by only based on
#' incidence difference (function \code{gci3}). \cr
#'
#' @inheritParams gci2
#' @param ptnum positive integer, number of AEs to be selected or plotted, default is 10
#' @param hierdata output from function \code{\link{Hier}}
#' @param isingdata output from function \code{\link{Ising}}
#' @param param a string, either "odds ratio" or "risk difference", indicate which summary statistic to be based on to plot the top AEs,
#' default is "risk difference"
#'
#' @seealso
#' \code{\link{preprocess}}, \code{\link{Hier}}, \code{\link{Ising}}
#'
#' @examples
#' \dontrun{
#' data(ADAE)
#' data(ADSL)
#' AEdata<-preprocess(adsl=ADSL, adae=ADAE)
#'
#' # run the Hierarchical model
#' INITS1<-list(mu.gamma.0=0.1, tau.gamma.0=0.1, mu.theta.0=0.1, tau.theta.0=0.1, alpha.pi=2, beta.pi=2)
#' INITS2<-list(mu.gamma.0=1, tau.gamma.0=1, mu.theta.0=1, tau.theta.0=1, alpha.pi=10, beta.pi=10)
#' INITS <- list(INITS1,INITS2)
#' HIERDATA<-Hier(aedata=AEdata, inits=INITS, n_burn=1000, n_iter=1000, thin=20, n_adapt=1000, n_chain=2)
#'
#' # run the Ising model
#' RHO<-rep(1,dim(AEdata)[1])
#' THETA<-0.02
#' SIM<-c(5000,1000,20)
#' BETA.AB<-c(0.25, 0.75)
#' ISINGDATA<-Ising(aedata = AEdata, beta.ab = BETA.AB, rho = RHO, theta = THETA, sim = SIM)
#'
#' gci3(aedata=AEdata)
#' HIPLOT(hierdata=HIERDATA, isingdata=ISINGDATA, aedata=AEdata)
#' HIPLOT(hierdata=HIERDATA, isingdata=ISINGDATA, aedata=AEdata, ptnum=15, param="odds ratio")
#' }
#'
#'
#'
#' @export
#'


gci3<-function(aedata, ptnum=10, conf.level=0.95){

  # This function calculate the incidence rate for each adverse event (AE) in treatment
  # and in control group and also the confidence interval for this incidence getting
  # from fisher exact test
  # then the PT of AEs with the highest ptnum incidence rate difference between treatment and control
  # group were returned

  # The input aedata is the output from preprocess,
  # it has the follwing fixed variable names
  # the data includes AE, not total subjects
  # Nc: number of control subjects
  # AEc: number of subjects with AE, for each term, out of Nc subjects, control group
  # Nt: number of treatment group subjects
  # AEt: number of subjects with AE, for each term, out of Nt subjects, treatment group
  # AEBODSYS: SOC of the adverse event
  # AEDECOD: PT of the adverse event
  # ptnum: the number of pt we want to return
  # conf.level: confident level used for fisher exact test with default be 0.95

  library(binom)
  aedata$mc = aedata$AEc/aedata$Nc
  aedata$llc = binom.confint(aedata$AEc,aedata$Nc, conf.level, method="exact")$lower
  aedata$ulc = binom.confint(aedata$AEc,aedata$Nc, conf.level, method="exact")$upper

  aedata$mt = aedata$AEt/aedata$Nt
  aedata$llt = binom.confint(aedata$AEt,aedata$Nt, conf.level, method="exact")$lower
  aedata$ult = binom.confint(aedata$AEt,aedata$Nt, conf.level, method="exact")$upper
  aedata$d = abs(aedata$mt - aedata$mc)
  aedata1=aedata[order(-aedata$d),]

  pdat= data.frame(lapply(aedata1, as.character), stringsAsFactors=FALSE)
  return(head(pdat, ptnum)$AEDECOD)
}


#########################################################################################################
#########################################################################################################

#' @rdname HierIsingplot
#'
#' @export
HIPLOT<-function(hierdata,isingdata, aedata, ptnum=10, param="risk difference" ){
  # hierdata is the result from function Hier
  # isingdata is the result from function Ising
  # aedata is the result from function preprocess
  # PTNUM is the number of AE we want to plot
  # param is the summary statistic we use to select the AE, it can be either "risk difference" or "odds ratio"
  # for this function, it will first select the top PTNUM number of AE based on the summary statistic param from the selected method
  # then it will plotted the mean, 2.5% quantile, 97.5% quantile of the param of these AE from both two models(methods)
  # it also has one more parameter aedata
  # this is for fucntion gci2 to get the top AEs with risk difference based on fisher exact test

  library(ggplot2)
  library(data.table)

  hierdata.new<-hierdata
  hierdata.new$gammaeq1<-NA
  inputdata<- rbind(hierdata.new,isingdata)

  # get the PTs of AE with top ptnum difference in risk difference based on fisher exact test
  FET<-gci3(aedata, ptnum)

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
    PT_Union<-intersect(PT_Hier, PT_Is)
    PT_Hier_Uni<-setdiff(PT_Hier, PT_Union)
    PT_Is_Uni<-setdiff(PT_Is, PT_Union)

    # get the number of PTs in each set
    n_Union<-length(PT_Union)
    n_Hier_Uni<-length(PT_Hier_Uni)
    n_Is_Uni<-length(PT_Is_Uni)

    # get the indicator for PT in different set
    # 0 for PT in union, 1 for PT unique in Hierarchical model
    # 2 for PT in Ising model
    test_Union<-inputdata[which(inputdata$PT %in% PT_Union),]
    test_Union$Ind<-"A"
    test_Union$Set<-"Intersect of both models"
    test_Hier_Uni<-inputdata[which(inputdata$PT %in% PT_Hier_Uni),]
    test_Hier_Uni$Ind<-"B"
    test_Hier_Uni$Set<-"Unique in Hierarchical Model"
    test_Is_Uni<-inputdata[which(inputdata$PT %in% PT_Is_Uni),]
    test_Is_Uni$Ind<-"C"
    test_Is_Uni$Set<-"Unique in Ising prior"

    # combine the dataset together
    test<-rbind(test_Union, test_Hier_Uni)
    test<-rbind(test, test_Is_Uni)

    # add one more column "category" for the indication of whether the AE is inside FET, AE slected by fisher exact test
    test$category<-0
    if(length(intersect(FET, test$PT)))  test[which(test$PT %in% FET), ]$category<-1
    # the if statement is to make sure there is at least one PT that in both FET and test

    # sort the test by Method, Ind, and -Diff_mean
    test<-test[order(test$Method, test$Ind, -test$Diff_mean),]
    # create the x column for plotting with Hierarchical model has x value be n-0.1
    # and Ising model has x value of n+0.1
    len<-dim(test)[1]/2
    PT_lable<-test$PT[1:len]
    X1<-seq(0.9, len-0.1, by=1)
    X2<-seq(1.1, len+0.1, by=1)
    X<-c(X1, X2)
    test$x<-X

    setnames(test, old=c("OR_2.5%", "OR_97.5%", "Diff_2.5%","Diff_97.5%" ),
             new=c("OR_L", "OR_U","Diff_L","Diff_U"))

    p <- ggplot(test, aes(x=x, y=Diff_mean, color=Set, shape=Method)) + geom_pointrange(aes(ymin=Diff_L, ymax=Diff_U))
    p1 <-p + labs(x = "Prefered Term",y="Risk difference",title = paste0("Top ", ptnum, " AE of mean risk difference"))

    a <- ifelse(test$category == 0, "red", "black")
    p2 <- p1 + theme(axis.text.x = element_text(color = a, size = 8, angle = 60, hjust = 1), axis.text.y = element_text(color = "black", size = 8))
    p3 <- p2 + scale_x_continuous(breaks=seq(1,len,by=1),labels=PT_lable)
    return(p3)
  }

  if (param=="odds ratio"){
    # first to get the top 10 AEs from each method
    test1_Hier<-subset(inputdata, Method=="Bayesian Hierarchical Model")
    test2_Hier<-test1_Hier[order(test1_Hier[, "OR_mean"], decreasing=TRUE), ]
    test3_Hier<-head(test2_Hier, ptnum)

    test1_Is<-subset(inputdata, Method=="Ising prior")
    test2_Is<-test1_Is[order(test1_Is[, "OR_mean"], decreasing=TRUE), ]
    test3_Is<-head(test2_Is, ptnum)

    # get the PT selected by each method
    PT_Hier<-test3_Hier$PT
    PT_Is<-test3_Is$PT

    # get the union of PTs, and PTs selected unique by each method
    PT_Union<-intersect(PT_Hier, PT_Is)
    PT_Hier_Uni<-setdiff(PT_Hier, PT_Union)
    PT_Is_Uni<-setdiff(PT_Is, PT_Union)

    # get the number of PTs in each set
    n_Union<-length(PT_Union)
    n_Hier_Uni<-length(PT_Hier_Uni)
    n_Is_Uni<-length(PT_Is_Uni)

    # get the indicator for PT in different set
    # 0 for PT in union, 1 for PT unique in Hierarchical model
    # 2 for PT in Ising model
    test_Union<-inputdata[which(inputdata$PT %in% PT_Union),]
    test_Union$Ind<-"A"
    test_Union$Set<-"Intersect of both models"
    test_Hier_Uni<-inputdata[which(inputdata$PT %in% PT_Hier_Uni),]
    test_Hier_Uni$Ind<-"B"
    test_Hier_Uni$Set<-"Unique in Hierarchical Model"
    test_Is_Uni<-inputdata[which(inputdata$PT %in% PT_Is_Uni),]
    test_Is_Uni$Ind<-"C"
    test_Is_Uni$Set<-"Unique in Ising prior"

    # combine the dataset together
    test<-rbind(test_Union, test_Hier_Uni)
    test<-rbind(test, test_Is_Uni)

    # add one more column "category" for the indication of whether the AE is inside FET, AE slected by fisher exact test
    test$category<-0
    if(length(intersect(FET, test$PT)))  test[which(test$PT %in% FET), ]$category<-1
    # the if statement is to make sure there is at least one PT that in both FET and test

    # sort the test by Method, Ind, and -Diff_mean
    test<-test[order(test$Method, test$Ind, -test$OR_mean),]

    # create the x column for plotting with Hierarchical model has x value be n-0.1
    # and Ising model has x value of n+0.1
    len<-dim(test)[1]/2
    PT_lable<-test$PT[1:len]
    X1<-seq(0.9, len-0.1, by=1)
    X2<-seq(1.1, len+0.1, by=1)
    X<-c(X1, X2)
    test$x<-X

    setnames(test, old=c("OR_2.5%", "OR_97.5%", "Diff_2.5%","Diff_97.5%" ),
             new=c("OR_L", "OR_U","Diff_L","Diff_U"))

    p <- ggplot(test, aes(x=x, y=OR_mean, color=Set, shape=Method)) + geom_pointrange(aes(ymin=OR_L, ymax=OR_U))
    p1 <-p + labs(x = "Prefered Term",y="Odds ratio",title = paste0("Top ", ptnum, " AE of mean odds ratio"))

    a <- ifelse(test$category == 0, "red", "black")
    p2 <- p1 + theme(axis.text.x = element_text(color = a, size = 8, angle = 60, hjust = 1), axis.text.y = element_text(color = "black", size = 8))
    p3 <- p2 + scale_x_continuous(breaks=seq(1,len,by=1),labels=PT_lable)
    return(p3)
  }
}
