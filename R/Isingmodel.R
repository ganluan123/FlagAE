#########################################################################################################
#########################################################################################################

#' @name Isingpriormodel
#'
#' @title Bayesian Model with Ising Latent Variables
#'
#' @description Functions here are to take the orginized data (output from
#'   \code{preprocess2}) and apply the Bayesian model with Ising latent variables.
#'   See details for model description and difference between each function.
#'
#' @param aedata output from function \code{\link{preprocess2}}
#' @param beta.ab numeric vector with length of 2, is the prior for beta distribution
#' @param rho numeric vector with length equals to the number of rows of data frame aedata
#' @param theta numeric, \code{rho} and \code{theta} are parameters for Ising prior
#' @param sim numeric vecotr with length of 3, integer for each element, sim[1] is the
#' number of iterations of brun in; sim[2] is the number of interactions to recorded;
#' sim[3] is like the parameter thin in MCMC settings. The total number of iterations running is
#' \eqn{sim[1]+sim[2]*sim[3]}
#' @param isingraw output from function \code{\link{Ising_history}}
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
#' then summarize the result for each AE and merge the result with the raw data.
#' \strong{\code{Isinggetpi}}: \cr
#' This function calculates pit (incidence of AE in treatment group) and
#' pic (incidence of AE in control group) from the output of \code{Ising_history}
#  The output is used for Loss function.
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
#' AEdata<-preprocess2(adsl=ADSL, adae=ADAE, TreatCol="TREATMENT", drug="xyz")
#' RHO<-rep(1,dim(AEdata)[1])
#' THETA<-0.02
#' SIM<-c(5000,1000,20)
#' BETA.AB<-c(0.25, 0.75)
#' ISINGRAW<-Ising_history(aedata = AEdata, beta.ab = BETA.AB, rho = RHO, theta = THETA, sim = SIM)
#' SUM_ISING<-sum_Ising(ISINGRAW)
#' ISINGMODEL<-Ising(aedata = AEdata, beta.ab = BETA.AB, rho = RHO, theta = THETA, sim = SIM)
#' ISINGPI<-Isinggetpi(aedata = AEdata, isingraw=ISINGRAW)
#' }
#'
#' @seealso
#' \code{\link{preprocess2}}
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

Ising_history <- function (aedata, beta.ab, rho, theta, sim){

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

  # parameter beta.ab is a two element vector, which is the prior for beta distribution
  # parameter rho, and theta are the hyperparameters for Ising prior

  # sim is a vector with 3 elements, the first element of sim is the iterations of burn in,
  # the second element of sim is the number of iterations we want to record
  # the third element of sim is the thin
  # please note the total number of iteration runing is sim[1]+sim[2]*sim[3]

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
  a.t<-beta.ab[1];b.t<-beta.ab[2];a.c<-beta.ab[1];b.c<-beta.ab[2]

  #simulation parameters
  burn.in<-sim[1]; T<-sim[2]; thin<-sim[3]

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
    logM0<-lbeta(a.t+y[k], N-y[k]+b.t)-lbeta(a.t, b.t)
    #model 1:treatment effect-->different pi across treatment arm
    #-->gamma.k=0
    logM1<-lbeta(a.t+y.t[k], n.t-y.t[k]+b.t)-lbeta(a.t, b.t)+
      lbeta(a.c+y.c[k], n.c-y.c[k]+b.c)-lbeta(a.c, b.c)
    logdiff<-logM0-logM1
    BF[k]<-exp(logdiff)
  }


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
        pi.t.m[t,k]<- rbeta(1, a.t+y.m[k], N-y.m[k]+b.t)*g.t.m[k] +
          rbeta(1, a.t+y.t.m[k], n.t-y.t.m[k]+b.t)*(1-g.t.m[k])
        #set.seed(1)
        pi.c.m[t,k]<- pi.t.m[t,k]*g.t.m[k] +
          rbeta(1, a.c+y.c.m[k], n.c-y.c.m[k]+b.t)*(1-g.t.m[k])
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
  log.or <- apply(log(isingraw$pi.t/(1-isingraw$pi.t))-log(isingraw$pi.c/(1-isingraw$pi.c)),2,mean) # log odds ratio
  log.or.LL <- apply(log(isingraw$pi.t/(1-isingraw$pi.t))-log(isingraw$pi.c/(1-isingraw$pi.c)),2,quantile,probs=0.025) # lower 2.5% quantile of log odds ratio
  log.or.UL <- apply(log(isingraw$pi.t/(1-isingraw$pi.t))-log(isingraw$pi.c/(1-isingraw$pi.c)),2,quantile,probs=0.975) # upper 2.5% quantile of log odds ratio
  rd <- apply(isingraw$pi.t-isingraw$pi.c,2,mean) # difference of pi.t and pi.c, incidence rate difference
  rd.LL <- apply(isingraw$pi.t-isingraw$pi.c,2,quantile,probs=0.025) # lower 2.5% quantile of incidence rate difference
  rd.UL <- apply(isingraw$pi.t-isingraw$pi.c,2,quantile,probs=0.975) # upper 2.5% quantile of incidicen rate difference

  cbind(p.gamma.k.eq1,pT.gt.pC,log.or,log.or.LL,log.or.UL,rd,rd.LL,rd.UL)
}


#########################################################################################################
#########################################################################################################

#' @rdname Isingpriormodel
#'
#' @export
Ising <- function (aedata, beta.ab, rho, theta, sim) {
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
  R_Ising_history<-Ising_history(aedata=aedata, beta.ab=beta.ab, rho=rho, theta=theta, sim=sim)

  # summarize the primary result with sum_Ising
  out <- round(sum_Ising(R_Ising_history),3)
  Ising.out <- cbind(as.character(aedata$AEBODSYS),as.character(aedata$AEDECOD),out)
  Ising.plot <- as.data.frame(Ising.out)
  Ising.plot$SoC <- Ising.plot$V1
  Ising.plot$PT <- Ising.plot$V2

  # get the probability that gamma=1, ie. no difference between treatment and control group
  Ising.plot$gammaeq1<-as.numeric(as.character(Ising.plot$p.gamma.k.eq1))

  # get Odds ratio mean and quantile
  Ising.plot$OR_mean <- exp(as.numeric(as.character(Ising.plot$log.or)))
  Ising.plot$'OR_2.5%' <- exp(as.numeric(as.character(Ising.plot$log.or.LL)))
  Ising.plot$'OR_97.5%' <- exp(as.numeric(as.character(Ising.plot$log.or.UL)))

  # get risk disfference mean and quantile
  Ising.plot$Diff_mean <- as.numeric(as.character(Ising.plot$rd))
  Ising.plot$'Diff_2.5%' <- as.numeric(as.character(Ising.plot$rd.LL))
  Ising.plot$'Diff_97.5%' <- as.numeric(as.character(Ising.plot$rd.UL))

  Ising.plot$Method <- 'Ising prior'

  Ising.plot <- Ising.plot[,c('SoC','PT','Diff_mean','Diff_2.5%','Diff_97.5%','OR_mean','OR_2.5%','OR_97.5%','Method','gammaeq1')]

  # get the raw data from aedata
  Raw<-aedata
  names(Raw)[1:2]<-c("SoC", "PT")
  # merge with raw data
  Ising.plot<-merge(Raw, Ising.plot, by=c("SoC", "PT"))
}


#########################################################################################################
#########################################################################################################

#' @rdname Isingpriormodel
#'
#' @export
Isinggetpi<-function (aedata, isingraw){
  # this function is to get the pit and pic from the output of Ising_history
  # the output is used for Loss function
  # aedata is the output of preprocess2
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

