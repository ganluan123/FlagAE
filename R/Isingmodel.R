#########################################################################################################
#########################################################################################################

#' @name Isingpriormodel
#'
#' @title Bayesian Model with Ising Latent Variables
#'
#' @description Functions here are to take the orginized data (output from
#'   \code{preprocess}) and apply the Bayesian model with Ising latent variables.
#'   See details for model description and difference between each function.
#'
#' @param aedata output from function \code{\link{preprocess}}
#' @param alpha_ numeric, is the prior for beta distribution, beta distribution for both treatment and control group, alpha parameter of beta distribution
#' @param beta_ numeric, is the prior for beta distribution, beta distribution for both treatment and control group, beta parameter of beta distribution
#' @param alpha.t numeric, is the prior for beta distribution, beta distribution for treatment group, alpha parameter of beta distribution
#' @param beta.t numeric, is the prior for beta distribution, beta distribution for  treatment group, beta parameter of beta distribution
#' @param alpha.c numeric, is the prior for beta distribution, beta distribution for control group, alpha parameter of beta distribution
#' @param beta.c numeric, is the prior for beta distribution, beta distribution for control group, beta parameter of beta distribution
#'
#' @param theta numeric, \code{rho} and \code{theta} are parameters for Ising prior
#' @param rho either a number or numeric vector with length equals to the number of rows of data frame aedata.
#' If it is a single number, then all adverse events use the same hyperparameter of rho. If it is a numeric vector, then each AE has
#' its own hyperparameter of rho, and the sequence of rho value for each AE should be the same as the sequence of AE in aedata (AE in
#' aedata should be ordered by b and j).
#' @param isingraw output from function \code{\link{Ising_history}}
#' @param isingdata output from function \code{\link{Ising}}
#' @param ptnum positive integer, number of AEs to be selected or plotted, default is 10
#' @param param a string, either "odds ratio" or "risk difference", indicate which summary statistic to be based on to plot the top AEs,
#' default is "risk difference"
#' @param n_burn number of burn in for Gibbs Sampling
#' @param n_iter number of interation for Gibbs Sampling
#' @param thin thin for Gibbs Samping, parameters are recorded every thin-th interation
#' @param OR_ylim a numeric vector of two elements, used to set y-axis limit for plotting based on "odds ratio"
#'
#' @details
#' \strong{Model}:\cr
#' The model is based on the paper McEvoy, Bradley W., Rajesh R. Nandy, and Ram C. Tiwari.
#' "Bayesian approach for clinical trial safety data using an Ising prior."
#' Biometrics 69.3 (2013): 661-672. \cr
#' \strong{\code{Ising_history}}:\cr
#' This function takes formatted Binomial data and
#' performs the Bayesian analysis with Ising prior.#'
#' \strong{\code{sum_Ising}}:\cr
#' This function takes the result from function \code{Ising_history} as input
#' and summarizes the result from \code{Ising_history}.\cr
#' \strong{\code{Ising}}:\cr
#' This function takes the same parameters as function \code{Ising_history}
#' this function will take first take the result from \code{Ising_history} and
#' then summarize the result for each AE and merge the result with the raw data,and also the Raw risk difference, Raw odds ratio calculated from raw data.
#' \strong{\code{Isinggetpi}}: \cr
#' This function calculates pit (incidence of AE in treatment group) and
#' pic (incidence of AE in control group) from the output of \code{Ising_history}
#  The output is used for Loss function.
#' \code{Isingplot} first selects the top \code{ptnum} (an integer) AE based on the selected statistic (either "odds ratio" or "risk difference").
#' Then it plots the mean, 2.5% quantile, 97.5% quantile of the selected statistic of these AE. It shows the PTs of these AE from same SOC in same
#' color.  \cr
#' \code{Isingtable} creates a table for the detailed information for AE plotted in \code{Isingplot}.
#'
#'
#' @return
#' \code{Ising_history} returns a list with 3 matries, gamma, pi.t, and pi.c.
#' Each row of these three matries correspond to one recoreded iteration.
#' And each column of these three matries correspond to each of the AE.\cr
#' \emph{gamma} matrix records the indicator of non-differential risk for each AE. \cr
#' \emph{pi.t} records the incidence risk of treatment group for each AE. \cr
#' \emph{pi.c} records the incidence risk of control group for each AE. \cr
#'
#' \code{sum_Ising} returns the summary statistics for each AE.(Note: all the apply function are apply
#' on column, thus will result the summary statistics for each AE). \cr
#'
#' \code{Ising} returns the summary statistics for each AE combining with the raw data. The summary statistics include
#' summary statistics for incidence rate difference (mean, 2.5\% and 97.5\% percentile),
#' summary statistics for odds ratio (mean, 2.5\% and 97.5\% percentile), and also the mean of gamma, with is the probability of gamma=1.
#' The other columns from raw data include SoC, PT, Nt, Nc, AEt, AEc. \cr
#'
#' \strong{\code{Isinggetpi}}: \cr
#' This function calculates pit (incidence of AE in treatment group) and
#' pic (incidence of AE in control group) from the output of \code{Ising_history}. \cr
#'
#' @examples
#' \dontrun{
#' data(ADAE)
#' data(ADSL)
#' AEdata<-preprocess(adsl=ADSL, adae=ADAE)
#' RHO<-rep(1,dim(AEdata)[1])
#' ISINGRAW<-Ising_history(aedata = AEdata, n_burn=1000, n_iter=5000, thin=20, rho=1, theta=0.02)
#' ISINGRAW2<-Ising_history(aedata = AEdata, n_burn=1000, n_iter=5000, thin=20, alpha_=0.5, beta_=0.5,
#'                            alpha.t=0.5, beta.t=0.5, alpha.c=0.25, beta.c=0.75, rho=RHO, theta=0.02)
#' SUM_ISING<-sum_Ising(ISINGRAW)
#' ISINGDATA<-Ising(aedata = AEdata, n_burn=1000, n_iter=5000, thin=20, rho=1, theta=0.02)
#' ISINGDATA<-Ising(aedata = AEdata, n_burn=1000, n_iter=5000, thin=20, alpha_=0.5, beta_=0.5,
#'                            alpha.t=0.5, beta.t=0.5, alpha.c=0.25, beta.c=0.75, rho=RHO, theta=0.02)
#' ISINGPI<-Isinggetpi(aedata = AEdata, isingraw=ISINGRAW)
#'
#' Isingplot(ISINGDATA)
#' Isingplot(ISINGDATA, ptnum=15, param="odds ratio", OR_ylim=c(1,10))
#' ISINGTABLE<-Isingtable(ISINGDATA)
#' ISINGTABLE2<-Isingtable(ISINGDATA, ptnum=15, param="odds ratio")
#'
#' Compareplot(ISINGDATA)
#' # user can use a very big number(bigger than total PTs in dataset) to plot out all the PTs
#' Compareplot(ISINGDATA, ptnum=5000, param='odds ratio')
#' }
#'
#' @seealso
#' \code{\link{preprocess}}
#'
#' @references
#' McEvoy, Bradley W., Rajesh R. Nandy, and Ram C. Tiwari.
#' "Bayesian approach for clinical trial safety data using an Ising prior."
#' Biometrics 69.3 (2013): 661-672.
#'
#' @note The incidence difference is calculated by (t-c) and the odds ratio is calculated by t(1-c)/c(1-t), where t,c are
#' the incidences of one AE for treatment and control group, respectively.\cr
#'
#' @export

Ising_history <- function (aedata, n_burn, n_iter, thin, alpha_=0.25, beta_=0.75, alpha.t=0.25, beta.t=0.75,
                           alpha.c=0.25, beta.c=0.75, rho, theta ){

  # this function is to perform the Bayesian analysis with Ising prior
  # it will output a list with 3 matries, gamma, pi.t, and pi.c
  # each row of these three matries correspond to one recoreded iteration.
  # and each column of these three matries correspond to each of the AE
  # gamma matrix records the indicator of non-differential risk for each AE
  # pi.t records the incidence risk of treatment group for each AE
  # pi.c records the incidence risk of control group for each AE

  # The parameter aedata has the follwing fixed variable names
  # Nc: number of control subjects
  # AEc: number of subjects with AE, for each term, out of Nc subjects, control group
  # Nt: number of treatment group subjects
  # AEt: number of subjects with AE, for each term, out of Nt subjects, treatment group
  # b: index for Soc

  # parameter beta.alpha, and beta.beta are the priors for beta distribution
  # parameter rho, and theta are the hyperparameters for Ising prior

  # n_burn, n_iter, and thin are Gibbs sampling parameters, see help document for details.

  n.SOC <- max(aedata$b)
  # SOClist is a list for the index of AE that belongs to the same SoC
  # the mth element of SOC is the vector of the row number of aedata with aedata$b == m,
  # ie SoC number == m
  # the way we create SOClist asks us to make sure aedata is ordered by b and j
  aedata<-aedata[order(aedata$b, aedata$j),]

  SOClist<-list()
  track<-0
  for (m in 1:n.SOC){
    SOCmlength<-dim(aedata[aedata$b==m,])[1]
    SOClist[[m]]<-(track+1):(track+SOCmlength)
    track<-track+SOCmlength
  }

  n.t<-aedata$Nt[1]
  n.c<-aedata$Nc[1]
  N<-n.t+n.c
  y.t<-c(aedata$AEt)
  y.c<-c(aedata$AEc)
  y<-y.t+y.c
  # y is the sum of y.t and y.c for each PT


  K<-length(y.t)
  #beta
  # a.t<-beta.ab[1];b.t<-beta.ab[2];a.c<-beta.ab[1];b.c<-beta.ab[2]

  #simulation parameters
  burn.in<-n_burn; T<-floor(n_iter/thin); thin<-thin

  BF<-0
  #############################################################
  # note to avoid the overflow and underflow problem,
  # remove the choose() function in both MO and M1,
  # since it will cancel out after dividsion
  # also we use lbeta() function instead of beta() function
  # to avoid the underflow problem
  #############################################################

  for (k in 1:K) {
    #model 0: no treatment effect--> same pi across treatment arms
    #-->gamma.k=1
    logM0<-lbeta(alpha_+y[k], N-y[k]+beta_)-lbeta(alpha_, beta_)
    #model 1:treatment effect-->different pi across treatment arm
    #-->gamma.k=0
    logM1<-lbeta(alpha.t+y.t[k], n.t-y.t[k]+beta.t)-lbeta(alpha.t, beta.t)+
      lbeta(alpha.c+y.c[k], n.c-y.c[k]+beta.c)-lbeta(alpha.c, beta.c)
    logdiff<-logM0-logM1
    BF[k]<-exp(logdiff)
  }

  # create the vector for rho if rho is a single number
  if(length(rho)==1) rho<-rep(rho, dim(aedata)[1])

  #setup parallel backend to use multiple processors
  library(foreach)
  library(doParallel)
  cores=detectCores()
  cl<-makeCluster(cores[1]-1)
  registerDoParallel(cl)

  List<-list()

  List<-foreach(m = 1:n.SOC) %dopar%{

    rownum.m<-SOClist[[m]] # this is indexs (i in aedata) for SOC (b in aedata) equals to m
    K.m<-length(rownum.m)
    gamma.m<-matrix(0,nrow=T, ncol=K.m)
    g.t.m<-rep(0,K.m)
    pi.t.m <-matrix(0,nrow=T,ncol=K.m)
    pi.c.m<-pi.t.m
    #pi.h.m<-matrix(0, nrow=T, ncol=K.m)
    BF.m<-BF[rownum.m]
    rho.m<-rho[rownum.m]
    y.t.m<-y.t[rownum.m]
    y.c.m<-y.c[rownum.m]
    y.m<-y[rownum.m]

    # first the burn process
    for (i in 1:burn.in){

      o.m<-rank(runif(K.m))
      # o.m<-1:K.m
      for (k in 1:K.m){
        h<-1/BF.m[o.m[k]]*exp(-rho.m[o.m[k]]+theta*sum(1-2*g.t.m) - theta*(1-2*g.t.m[o.m[k]])) # excluding itself
        #set.seed(1)
        g.t.m[o.m[k]]<-rbinom(1,1,1/(1+h))
      }
    }


    # then the real running phase
    for (t in 1:T) {
      for (tt in 1:thin) {

        o.m<-rank(runif(K.m))
        # o.m<-1:K.m
        for (k in 1:K.m) {
          h<-1/BF.m[o.m[k]]*exp(-rho.m[o.m[k]]+theta*sum(1-2*g.t.m) - theta*(1-2*g.t.m[o.m[k]]))
          #set.seed(1)
          g.t.m[o.m[k]]<-rbinom(1,1,1/(1+h))
          #pi.h.m[t,k]<-1/(1+h)
        }
      }

      for (k in 1:K.m) {
        #set.seed(1)
        pi.t.m[t,k]<- rbeta(1, alpha_+y.m[k], N-y.m[k]+beta_)*g.t.m[k] +
          rbeta(1, alpha.t+y.t.m[k], n.t-y.t.m[k]+beta.t)*(1-g.t.m[k])
        #set.seed(1)
        pi.c.m[t,k]<- pi.t.m[t,k]*g.t.m[k] +
          rbeta(1, alpha.c+y.c.m[k], n.c-y.c.m[k]+beta.c)*(1-g.t.m[k])
      }

      gamma.m[t,]<-g.t.m
    }
    # List[[m]]<-list(gamma=gamma.m, pi.t=pi.t.m, pi.c=pi.c.m)
    List[[m]]<-list(gamma=gamma.m, pi.t=pi.t.m, pi.c=pi.c.m)

  }
  # stop the setting for parallel computing
  stopCluster(cl)


  # combine the result from each SOC to generate the final result
  gamma=List[[1]]$gamma
  pi.t=List[[1]]$pi.t
  pi.c=List[[1]]$pi.c
  #pi.h=List[[1]]$pi.h


  for (l in 2:n.SOC){
    gamma<-cbind(gamma, List[[l]]$gamma)
    pi.t<-cbind(pi.t, List[[l]]$pi.t)
    pi.c<-cbind(pi.c, List[[l]]$pi.c)
    #pi.h<-cbind(pi.h, List[[l]]$pi.h)

  }

  #Finallist<-list(gamma=gamma, pi.t=pi.t, pi.c=pi.c)
  Finallist<-list(gamma=gamma, pi.t=pi.t, pi.c=pi.c)
}


#########################################################################################################
#########################################################################################################

#' @rdname Isingpriormodel
#'
#' @export
sum_Ising <- function(isingraw) {
  # this function takes the result from function Ising_history as input
  # this function summarize the result from Ising_history
  # it will return the summary statistics for each AE
  # all the apply function are apply on column, thus will result the summary statistics for each AE
  p.gamma.k.eq1 <-apply(isingraw$gamma,2,mean) # probability that gamma=1
  pT.gt.pC <-apply(isingraw$pi.t>isingraw$pi.c,2,mean) # probability that pi.t > pi.c
  #log.rr <- apply(log(isingraw$pi.t)-log(isingraw$pi.c),2,mean) # difference of log(pi.t) and log(pi.c)
  log.or <- apply(log(isingraw$pi.t/(1-isingraw$pi.t))-log(isingraw$pi.c/(1-isingraw$pi.c)),2,quantile,probs=0.5) # log odds ratio
  log.or.LL <- apply(log(isingraw$pi.t/(1-isingraw$pi.t))-log(isingraw$pi.c/(1-isingraw$pi.c)),2,quantile,probs=0.025) # lower 2.5% quantile of log odds ratio
  log.or.UL <- apply(log(isingraw$pi.t/(1-isingraw$pi.t))-log(isingraw$pi.c/(1-isingraw$pi.c)),2,quantile,probs=0.975) # upper 2.5% quantile of log odds ratio
  rd <- apply(isingraw$pi.t-isingraw$pi.c,2,mean) # difference of pi.t and pi.c, incidence rate difference
  # rd <- apply(isingraw$pi.t,2,mean)
  # rd <- apply(isingraw$pi.c,2,mean)
  # rd <- apply(isingraw$pi.t-isingraw$pi.c,2,quantile,probs=0.5)
  rd.LL <- apply(isingraw$pi.t-isingraw$pi.c,2,quantile,probs=0.025) # lower 2.5% quantile of incidence rate difference
  rd.UL <- apply(isingraw$pi.t-isingraw$pi.c,2,quantile,probs=0.975) # upper 2.5% quantile of incidicen rate difference

  cbind(p.gamma.k.eq1,pT.gt.pC,log.or,log.or.LL,log.or.UL,rd,rd.LL,rd.UL)
}


#########################################################################################################
#########################################################################################################

#' @rdname Isingpriormodel
#'
#' @export
Ising <- function (aedata, n_burn, n_iter, thin, alpha_=0.25, beta_=0.75, alpha.t=0.25, beta.t=0.75,
                   alpha.c=0.25, beta.c=0.75, rho, theta) {
  # this function takes the same parameters as function Ising_history
  # this function will take first take the result from Ising_history and
  # then summarize the result for each AE and merge the result with the raw data
  # the summary statistics including:
  # summary statistics for incidence rate difference (mean, 2.5% and 97.5% percentile)
  # summary statistics for odds ratio (mean, 2.5% and 97.5% percentile)
  # also the mean of gamma, with is the probability of gamma=1
  # the other columns including:
  # SoC, PT, Nt, Nc, AEt, AEc

  aedata<-aedata[order(aedata$b),]
  R_Ising_history<-Ising_history(aedata=aedata, n_burn=n_burn, n_iter=n_iter, thin=thin,
                                 alpha_=alpha_, beta_=beta_, alpha.t=alpha.t, beta.t=beta.t,
                                 alpha.c=alpha.c, beta.c=beta.c, rho=rho, theta=theta)

  # summarize the primary result with sum_Ising
  out <- round(sum_Ising(R_Ising_history),3)
  Ising.out <- cbind(as.character(aedata$AEBODSYS),as.character(aedata$AEDECOD),out)
  Ising.plot <- as.data.frame(Ising.out)
  Ising.plot$SoC <- Ising.plot$V1
  Ising.plot$PT <- Ising.plot$V2

  # get the probability that gamma=1, ie. no difference between treatment and control group
  Ising.plot$gammaeq1<-as.numeric(as.character(Ising.plot$p.gamma.k.eq1))

  # get Odds ratio median and quantile
  Ising.plot$OR_median <- exp(as.numeric(as.character(Ising.plot$log.or)))
  Ising.plot$'OR_2.5%' <- exp(as.numeric(as.character(Ising.plot$log.or.LL)))
  Ising.plot$'OR_97.5%' <- exp(as.numeric(as.character(Ising.plot$log.or.UL)))

  # get risk disfference mean and quantile
  Ising.plot$Diff_mean <- as.numeric(as.character(Ising.plot$rd))
  Ising.plot$'Diff_2.5%' <- as.numeric(as.character(Ising.plot$rd.LL))
  Ising.plot$'Diff_97.5%' <- as.numeric(as.character(Ising.plot$rd.UL))

  Ising.plot$Method <- 'Ising prior'

  Ising.plot <- Ising.plot[,c('SoC','PT','Diff_mean','Diff_2.5%','Diff_97.5%','OR_median','OR_2.5%','OR_97.5%','Method','gammaeq1')]

  # get the raw data from aedata
  Raw<-aedata

  names(Raw)[1:2]<-c("SoC", "PT")
  # merge with raw data
  Ising.plot<-merge(Raw, Ising.plot, by=c("SoC", "PT"))
  Ising.plot <- Ising.plot[order(Ising.plot$SoC),]
  #remove the column 'b', 'i', 'j'
  drops<-c('b', 'i', 'j')
  Ising.plot<-Ising.plot[, !names(Ising.plot) %in% drops]
  # reorder the column
  Ising.plot[, c(1:7, 9:11, 8, 12:16)]

}


#########################################################################################################
#########################################################################################################

#' @rdname Isingpriormodel
#'
#' @export
Isinggetpi<-function (aedata, isingraw){
  # this function is to get the pit and pic from the output of Ising_history
  # the output is used for Loss function
  # aedata is the output of preprocess
  # isingraw is the output of Ising_history
  # it will output a list with two elements, pit and pic
  # pit contains the columns of SoC, PT for each AE and also the incidence rate
  # for AE in treatment group for each iteraction
  # pic contains the columns of SoC, PT for each AE and also the incidence rate
  # for AE in control group for each iteraction

  # the result from isingraw is ordered by b and j
  # we need to make sure aedata is also ordered by b and j
  aedata<-aedata[order(aedata$b, aedata$j),]
  # combine SoC and PT together with incidence rate
  pit<-cbind(as.character(aedata$AEBODSYS), as.character(aedata$AEDECOD), t(isingraw$pi.t))
  pic<-cbind(as.character(aedata$AEBODSYS), as.character(aedata$AEDECOD), t(isingraw$pi.c))
  # convert to data frame
  pit<-as.data.frame(pit)
  pic<-as.data.frame(pic)

  # rename of first two columns of pit and pic
  names(pit)[1:2]<-c("SoC", "PT")
  names(pic)[1:2]<-c("SoC", "PT")

  List<-list(pit=pit, pic=pic)
  return (List)
}

#########################################################################################################
#########################################################################################################

#' @rdname Isingpriormodel
#'
#' @export
Isingplot<-function(isingdata, ptnum=10, param="risk difference", OR_xlim=c(0,5) ){
  # isingdata is the result from function Ising
  # ptnum is the number of AE we want to plot
  # param is the summary statistic we use to select the AE, it can be either "risk difference" or "odds ratio"
  # for this function, it will first select the top ptnum number of AE based on the summary statistic param from the selected method
  # then it will plotted the mean, 2.5% quantile, 97.5% quantile of the param of these AE


  library(ggplot2)
  library(data.table)
  inputdata<-copy(isingdata)

  if (param=="risk difference"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "Diff_mean"], decreasing=TRUE), ]
    ptnum<-min(ptnum, dim(test)[1])
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
    p<-p+geom_text(aes(label=textAE, x=New_Diff+0.1, y=yloc), size=3, show.legend = FALSE, na.rm = TRUE)

    # add title
    p<-p+ggtitle(paste0("Top ", ptnum, " AE of mean risk difference \nplotted with 95% credible interval"))
    p<-p + theme(plot.title = element_text(size=15, hjust=0.5))

    # ylable and x label
    p<-p+xlab("Risk Difference")+ylab("PT")
    p<-p + theme(axis.title.x = element_text(size=13)) + theme(axis.title.y = element_text(size=13))

    # x-axis coordinate and y-axis coordinate
    p<-p + theme(axis.text.y = element_text(color = as.factor(test$SoC[1:ptnum]), size = 13,angle=60, hjust = 1))
    p<-p + scale_y_continuous(breaks=seq(from=ptnum,to=1,by=-1),labels=test$PT[1:ptnum])
    p<-p + theme(axis.text.x=element_text(size=13))

    # legend size
    p<-p+theme(legend.text = element_text(size=13), legend.title = element_blank(), legend.position = "bottom"
               , legend.box = 'vertical', legend.background = element_rect(fill="lightblue"))

    return(p)

  }

  if (param=="odds ratio"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "OR_median"], decreasing=TRUE), ]
    ptnum<-min(ptnum, dim(test)[1])
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
    p<-p+ggtitle(paste0("Top ", ptnum, " AE of median odds ratio \nplotted with 95% credible interval"))
    p<-p + theme(plot.title = element_text(size=15, hjust=0.5))

    # ylable and x label
    p<-p+xlab("Odds ratio")+ylab("PT")
    p<-p + theme(axis.title.x = element_text(size=13)) + theme(axis.title.y = element_text(size=13))

    # x-axis coordinate and y-axis coordinate
    p<-p + theme(axis.text.y = element_text(color = as.factor(test$SoC[1:ptnum]), size = 13,angle=60, hjust = 1))
    p<-p + scale_y_continuous(breaks=seq(from=ptnum,to=1,by=-1),labels=test$PT[1:ptnum])
    p<-p + theme(axis.text.x=element_text(size=13))

    # legend size
    p<-p+theme(legend.text = element_text(size=13), legend.title = element_blank(), legend.position = "bottom"
               , legend.box = 'vertical', legend.background = element_rect(fill="lightblue"))

    return(p)
  }
}



#########################################################################################################
#########################################################################################################

#' @rdname Isingpriormodel
#'
#' @export
Isingtable<-function(isingdata, ptnum=10, param="risk difference" ){
  # this function takes the same parameter as function Hierplot
  # it returns the detailed information of AEs plotted in Hierplot

  library(data.table)
  inputdata<-isingdata

  if (param=="risk difference"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "Diff_mean"], decreasing=TRUE), ]
    ptnum<-min(ptnum, dim(test)[1])
    test<-head(test, ptnum)
  }

  if (param=="odds ratio"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "OR_median"], decreasing=TRUE), ]
    ptnum<-min(ptnum, dim(test)[1])
    test<-head(test, ptnum)
  }
  return(test)
}

