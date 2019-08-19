#########################################################################################################
#########################################################################################################

#' @name Hiermodel
#'
#' @title Three stages Bayesian Hierarchical Model
#'
#' @description  Functions here are to take the orginized data (output from
#'   \code{preprocess}) and apply the three stages Bayesian Hierarchical Model.
#'   See details for model description and difference between each function.
#'
#' @param aedata output from function \code{\link{preprocess}}
#' @param n_burn integer, number of interations without saving posterior samples
#' @param n_iter integer, number of interations saving posterior samples with every \code{thin}th iteration
#' @param thin integer, samples are saved for every \code{thin}th iteration
#' @param n_adapt integer, number of adaptations
#' @param n_chain number of MCMC chains
#' @param alpha.gamma... alpha.gamma and other parameters are hyperparameters for prior distribution,
#' see the reference for the meaning of each parameters
#' @param hierdata utput from function \code{\link{Hier}}
#' @param ptnum positive integer, number of AEs to be selected or plotted, default is 10
#' @param param a string, either "odds ratio" or "risk difference", indicate which summary statistic to be based on to plot the top AEs,
#' default is "risk difference"
#' @param OR_ylim a numeric vector of two elements, used to set y-axis limit for plotting based on "odds ratio"
#'
#'
#' @details \strong{Model}: \cr Here the 3-stage hierarchical bayesian model was
#'   used to model the probability of AEs. It is model 1b (Bayesian Logistic
#'   Regression Model with Mixture Prior on Log-OR) in H. Amy Xia , Haijun Ma &
#'   Bradley P. Carlin (2011) Bayesian Hierarchical Modeling for Detecting
#'   SaBCIy Signals in Clinical Trials, Journal of Biopharmaceutical Statistics,
#'   21:5, 1006-1029, DOI: 10.1080/10543406.2010.520181) \cr
#'   \strong{\code{Hier_history}}: \cr
#'   This function takes formatted Binomial data and
#'   output Gibbs sample of the defined parametes. The output is a dataframe
#'   with each column represent one parameter and each row is the output from
#'   one sampling/one iteration. Diff, OR, gamma, and theta are the parameters recorded.\cr
#'   \emph{Diff}: is the difference of incidence
#'   of AE between treatment and control group (treatment - control) \cr
#'   \emph{OR}: is the odds ratio for the incidence of AE between treatment and
#'   control group (treatment over control): t(1-c)/c(1-t), where t,c are the incidences of one AE for treatment
#'   and control group, respectively.\cr
#'   \emph{gamma}: logit(incidence of AE in control group) = gamma \cr
#'   \emph{theta}: logit(incidence of AE in treatment group) = gamma + theta; and OR = exp(theta) \cr
#'   The result for Diff, OR, gamma, and theta are ordered by j and then by b.
#'   For example the result is like Diff.1.1, Diff.2.1, Diff.3.1, Diff.4.1 Diff.1.2, Diff.2.2,
#'   Diff.3.2 and so on \cr
#'   \strong{\code{sum_Hier}}:\cr
#'   This function takes the output from \code{Hier_history} and return the summary
#'   statistics for each parameter recorded by \code{Hier_history}. The summary function is applied on each column. \cr
#'   \strong{\code{Hier}}:\cr
#'   This function takes the same input as \code{Hier_history} and calculate the summary statistics
#'   for output from \code{Hier_history}.
#'   It outputs the summary statistics for each AE, combining with raw data, and also the Raw risk difference, Raw odds ratio
#'   calculated from raw data.\cr
#'   \strong{\code{Hiergetpi}}: \cr
#'   This function calculates pit (incidence of AE in treatment group) and
#'   pic (incidence of AE in control group) from the output of \code{Hier_history}
#'   The output is used for Loss function. \cr
#'   \code{Hierplot} first selects the top \code{ptnum} (an integer) AE based on the selected statistic (either "odds ratio" or "risk difference").
#'   Then it plots the mean, 2.5% quantile, 97.5% quantile of the selected statistic of these AE. It shows the PTs of these AE from same SOC in same
#'   color.  \cr
#'   \code{Hiertable} creates a table for the detailed information for AE plotted in \code{Hierplot}.
#'
#'
#'
#' @return
#' \strong{\code{Hier_history}}\cr
#' It returns a dataframe with each column represent one parameter
#' and each row is the output from one sampling/one iteraction (like the output of \code{\link[R2jags]{coda.samples}}) \cr
#' \strong{\code{sum_Hier}} \cr
#' It returns the summary statistics for each parameter recorded by \code{Hier_history}. \cr
#' \strong{\code{Hier}} \cr
#' It returns the summary statistics for each AE, combining with raw data.
#' The summary statistics including:
#' summary statistics for incidence rate difference (mean, 2.5\% and 97.5\% percentile);
#' summary statistics for odds ratio (mean, 2.5\% and 97.5\% percentile).
#' The other columns include SoC, PT, Nt, Nc, AEt, and AEc.
#' \strong{\code{Hiergetpi}}: \cr
#' This function calculates pit (incidence of AE in treatment group) and
#' pic (incidence of AE in control group) from the output of \code{Hier_history}.\cr
#'
#' @examples
#' \dontrun{
#' data(ADAE)
#' data(ADSL)
#' AEdata<-preprocess(adsl=ADSL, adae=ADAE)
#' HIERRAW<-Hier_history(aedata=AEdata, n_burn=1000, n_iter=1000, thin=20, n_adapt=1000, n_chain=2)
#' HIERRAW2<-Hier_history(aedata=AEdata, n_burn=1000, n_iter=1000, thin=20, n_adapt=1000, n_chain=1)
#' # you can specify the hyperparameter as shown below
#' HIERRAW3<-Hier_history(aedata=AEdata, n_burn=1000, n_iter=1000, thin=20, n_adapt=1000, n_chain=1,
#' alpha.gamma=5, beta.gamma=1, alpha.theta=3, beta.theta=1, mu.gamma.0.0=0.1, tau.gamma.0.0=0.1, alpha.gamma.0.0=5,
#' beta.gamma.0.0=1, lambda.alpha=0.1, lambda.beta=0.1, mu.theta.0.0=0.1, tau.theta.0.0=0.1,alpha.theta.0.0=3, beta.theta.0.0=1)
#'
#' HIERDATA<-Hier(aedata=AEdata, n_burn=1000, n_iter=1000, thin=20, n_adapt=1000, n_chain=2)
#' HIERPI<-Hiergetpi(aedata=AEdata, hierraw=HIERRAW)
#'
#' Hierplot(HIERDATA)
#' Hierplot(HIERDATA, ptnum=15, param="odds ratio")
#' HIERTABLE<-Hiertable(HIERDATA)
#' HIERTABLE2<-Hiertable(HIERDATA, ptnum=15, param="odds ratio")
#' }
#'
#' @seealso
#' \code{\link{preprocess}}
#'
#' @references H. Amy Xia , Haijun Ma & Bradley P. Carlin (2011) Bayesian
#'   Hierarchical Modeling for Detecting Safety Signals in Clinical Trials,
#'   Journal of Biopharmaceutical Statistics, 21:5, 1006-1029, DOI:
#'   10.1080/10543406.2010.520181)
#'
#' @export

Hier_history<- function(aedata, n_burn, n_iter, thin, n_adapt, n_chain, alpha.gamma=3, beta.gamma=1,
                        alpha.theta=3, beta.theta=1, mu.gamma.0.0=0, tau.gamma.0.0=0.1, alpha.gamma.0.0=3,
                        beta.gamma.0.0=1, lambda.alpha=0.1, lambda.beta=0.1, mu.theta.0.0=0, tau.theta.0.0=0.1,
                        alpha.theta.0.0=3, beta.theta.0.0=1) {

  # This function takes formatted Binomial data and output
  # Gibbs sample of the defined parameters
  # the output is a dataframe with each column represent one parameter
  # each row is the output from one sampling/one iteration
  # the result for Diff, OR, gamma, and theta are ordered by j and then by b
  # for example the result is like Diff.1.1, Diff.2.1, Diff.3.1, Diff.4.1
  # Diff.1.2, Diff.2.2, Diff.3.2 and so on
  # aedata is the output from function preprocess
  # each row of the aedata corresponds to one AE
  # it contains the following columns:
  # Nc: number of total patients in control group
  # Nt: number of total patients in treatment group
  # b: SOC # of the AE
  # j: PT # of the AE
  # AEt: number of patients with AE in treatment group
  # AEc: number of patients with AE in control group
  ##################################################
  # n_burn, n_iter, thin, n_adapt, n_chain are for MCMC
  # n_burn: number of interations without saving posterior samples
  # n_iter: number of interations saving posterior samples with every nth iteration
  # n here is given by thin
  # n_adapt: number of adaptations
  # n_chain: number of MCMC chains


  #############################################
  ## M1b: Binomial model with mixture prior ###
  #############################################

  # create the model

  model.binom<-"model{
  for (i in 1:Nae) {
  X[i] ~ dbin(c[b[i], j[i]], Nc)
  Y[i] ~ dbin(t[b[i], j[i]], Nt)

  logit(c[b[i], j[i]]) <- gamma[b[i], j[i]]
  logit(t[b[i], j[i]]) <- gamma[b[i], j[i]] + theta[b[i], j[i]]

  gamma[b[i], j[i]] ~ dnorm(mu.gamma[b[i]], tau.gamma[b[i]])
  p0[i] ~ dbern(pi[b[i]] ) # prob of point mass
  theta1[b[i], j[i]] ~ dnorm(mu.theta[b[i]], tau.theta[b[i]])

  # theta=0 w.p. pi[i] and theta=theta1 w.p. 1-pi[i]

  theta[b[i], j[i]] <- (1- p0[i]) * theta1[b[i], j[i]]

  OR[b[i],j[i]] <- exp(theta[b[i],j[i]] )
  # ORpv2[b[i], j[i]] <- step(OR[b[i],j[i]] - 2 )  # OR >= 2

  # ORpv2[b[i], j[i]] <- step(OR[b[i],j[i]] -1.2 ) # OR >= 1.2

  # ORpv[b[i], j[i]] <- 1- step(-OR[b[i],j[i]]) # OR >1

  # RD[b[i], j[i]] <- t[b[i], j[i]] - c[b[i], j[i]]
  # RDpv[b[i], j[i]] <- 1 - step(c[b[i], j[i]] - t[b[i], j[i]])
  # RD>0
  # RDpv2[b[i], j[i]] <- step(t[b[i], j[i]] - c[b[i], j[i]]- 0.02) # RD>=2%
  # RDpv5[b[i], j[i]] <- step(t[b[i], j[i]] - c[b[i], j[i]]- 0.05) # RD>=5%

  D[i] <- X[i]*log(c[b[i], j[i]]) + (Nc-X[i])*log(1-c[b[i], j[i]]) + Y[i]*log(t[b[i], j[i]]) + (Nt-Y[i])*log(1-t[b[i],j[i]])


  Diff[b[i], j[i]] <- t[b[i], j[i]] - c[b[i], j[i]]
  }

  Dbar <- -2* sum(D[]) # -2logL without normalizing constant
  # SOC level parameters

  for(k in 1:B){
  pi[k] ~ dbeta(alpha.pi, beta.pi)
  mu.gamma[k] ~ dnorm(mu.gamma.0, tau.gamma.0)
  tau.gamma[k] ~ dgamma(alpha.gamma,beta.gamma)
  mu.theta[k] ~ dnorm(mu.theta.0, tau.theta.0)
  tau.theta[k] ~ dgamma(alpha.theta,beta.theta)
  }

  # hyperpriors for gammas;
  mu.gamma.0 ~ dnorm(mu.gamma.0.0, tau.gamma.0.0)
  tau.gamma.0 ~ dgamma(alpha.gamma.0.0, beta.gamma.0.0)

  # hyperpriors for thetas;
  mu.theta.0 ~ dnorm(mu.theta.0.0, tau.theta.0.0)
  tau.theta.0 ~ dgamma(alpha.theta.0.0,beta.theta.0.0)

  # hyperpriors for pi?s;
  alpha.pi ~ dexp(lambda.alpha)I(1,)
  beta.pi ~ dexp(lambda.beta)I(1,)
}"


  param<-c("OR", "Diff", "gamma", "theta")

  data <- list(Nae = nrow(aedata), Nc = aedata$Nc[1], Nt = aedata$Nt[1], B = max(aedata$b),
               b = aedata$b, j = aedata$j, Y = aedata$AEt, X = aedata$AEc,
               alpha.gamma=alpha.gamma, beta.gamma=beta.gamma,
               alpha.theta=alpha.theta, beta.theta=beta.theta, mu.gamma.0.0=mu.gamma.0.0, tau.gamma.0.0=tau.gamma.0.0,
               alpha.gamma.0.0=alpha.gamma.0.0,
               beta.gamma.0.0=beta.gamma.0.0, lambda.alpha=lambda.alpha, lambda.beta=lambda.beta,
               mu.theta.0.0=mu.theta.0.0, tau.theta.0.0=tau.theta.0.0,
               alpha.theta.0.0=alpha.theta.0.0, beta.theta.0.0=beta.theta.0.0)


  # we use parallel computing for n_chain>1
  if (n_chain>1){
    #setup parallel backend to use multiple processors
    library(foreach)
    library(doParallel)
    cores<-detectCores()
    cl<-makeCluster(cores[1]-1)
    registerDoParallel(cl)

    param.est<-list()
    param.est<-foreach(m=1:n_chain) %dopar% {
      library(mcmcplots)
      library(rjags)
      library(R2jags)
      temp.fit <- jags.model(textConnection(model.binom), data=data, n.chains=1, n.adapt=n_adapt, quiet=TRUE)
      update(temp.fit, n.iter=n_burn)

      # summary of posterior samples
      temp.param.samples <- coda.samples(temp.fit,param,n.iter=n_iter,thin=thin)
      temp.param.est <- data.frame(as.matrix(temp.param.samples))
      param.est[[m]]<-temp.param.est
    }
    stopCluster(cl)

    ## combine the result from seperate chains together
    Final.est<-param.est[[1]]
    for (m in 2:n_chain){
      Final.est<-rbind(Final.est, param.est[[m]])
    }
  }

  else {
    library(mcmcplots)
    library(rjags)
    library(R2jags)
    temp.fit <- jags.model(textConnection(model.binom),data=data,n.chains=1, n.adapt=n_adapt, quiet=TRUE)
    update(temp.fit, n.iter=n_burn)

    # summary of posterior samples
    temp.param.samples <- coda.samples(temp.fit,param,n.iter=n_iter,thin=thin)
    temp.param.est <- data.frame(as.matrix(temp.param.samples))
    Final.est<-temp.param.est
  }

  return(Final.est)
  }





#########################################################################################################
#########################################################################################################

#' @rdname Hiermodel
#' @export
sum_Hier <- function(hierraw){
  # this functin will take the output from Hier_history as input
  # and return the summary statistics for each parameter
  # the summary function is applyed on each column
  xbar <- mean(hierraw)
  xsd  <- sd(hierraw)
  x2.5 <- quantile(hierraw,0.025)
  x25  <- quantile(hierraw,0.25)
  xmdn <- quantile(hierraw,0.5)
  x75  <- quantile(hierraw,0.75)
  x97.5<- quantile(hierraw,0.975)
  out <- c(xbar,xsd,x2.5,x25,xmdn,x75,x97.5)
  return(out)
}


#########################################################################################################
#########################################################################################################

#' @rdname Hiermodel
#' @export
Hier<- function(aedata, n_burn, n_iter, thin, n_adapt, n_chain, alpha.gamma=3, beta.gamma=1,
                alpha.theta=3, beta.theta=1, mu.gamma.0.0=0, tau.gamma.0.0=0.1, alpha.gamma.0.0=3,
                beta.gamma.0.0=1, lambda.alpha=0.1, lambda.beta=0.1, mu.theta.0.0=0, tau.theta.0.0=0.1,
                alpha.theta.0.0=3, beta.theta.0.0=1){

  # this function take the same input as Hier_history
  # this function will get the summary statistics for output from Hier_history
  # it will give the summary statistics for each AE and combine with raw data
  # the summary statistics including:
  # summary statistics for incidence rate difference (mean, 2.5% and 97.5% percentile)
  # summary statistics for odds ratio (mean, 2.5% and 97.5% percentile)
  # the other columns including:
  # SoC, PT, Nt, Nc, AEt, AEc

  oest<-Hier_history(aedata=aedata, n_burn=n_burn, n_iter=n_iter, thin=thin, n_adapt=n_adapt, n_chain=n_chain, alpha.gamma=alpha.gamma, beta.gamma=beta.gamma,
                     alpha.theta=alpha.theta, beta.theta=beta.theta, mu.gamma.0.0=mu.gamma.0.0, tau.gamma.0.0=tau.gamma.0.0,
                     alpha.gamma.0.0=alpha.gamma.0.0, beta.gamma.0.0=beta.gamma.0.0, lambda.alpha=lambda.alpha,
                     lambda.beta=lambda.beta, mu.theta.0.0=mu.theta.0.0, tau.theta.0.0=tau.theta.0.0,
                     alpha.theta.0.0=alpha.theta.0.0, beta.theta.0.0=beta.theta.0.0)

  # Get the mean, standard devision, quantile of 2.5%, 25%, 50%, 75%, and 97.5% for the parameters interested
  # get the summary for parameter Dbar, Diff, and OR
  oest_Diff<-oest[, grepl("Diff", names(oest))]
  oest_OR<-oest[, grepl("OR", names(oest))]
  oest_DiffOR<-cbind(oest_Diff, oest_OR)

  param.sum <- sapply(oest_DiffOR,sum_Hier)
  summary <- as.data.frame(t(param.sum))

  # Associate the parameters with SOC and PT names
  # sort aedata by j, since parameters in summary are sorted by j and then by b.
  AEDECOD <- aedata[order(aedata$j, aedata$b),]
  SoC <- rep(as.character(AEDECOD$AEBODSYS),2)
  PT <- rep(as.character(AEDECOD$AEDECOD),2)
  Sub <- c(rep('Diff',nrow(aedata)),rep('OR',nrow(aedata)))
  col.name <- c("Mean","SD","2.5%","25%","50%","75%","97.5%")
  colnames(summary) <- col.name
  summary <- cbind(Sub,SoC,PT,summary)

  # Summary statistics for parameter Diff
  summary.diff <- summary[summary$Sub=='Diff',c('SoC','PT','Mean','2.5%','97.5%')]
  colnames(summary.diff)[3:5] <- c('Diff_mean','Diff_2.5%','Diff_97.5%')

  # Summary statistics for parameter OR
  summary.OR <- summary[summary$Sub=='OR',c('SoC','PT','50%','2.5%','97.5%')]
  colnames(summary.OR)[3:5] <- c('OR_median','OR_2.5%','OR_97.5%')

  # merge summary statistics with raw data
  out <- merge(summary.diff,summary.OR,by=c('SoC','PT')) # this merge function also sort the resulting dataframe by Soc and PT
  # get the raw data from aedata
  Raw<-aedata

  names(Raw)[1:2]<-c("SoC", "PT")
  out<-merge(Raw, out, by=c("SoC", "PT"))
  Hier.plot <- out[order(out$SoC),]
  Hier.plot$Method = 'Bayesian Hierarchical Model'

  #remove the column 'b', 'i', 'j'
  drops<-c('b', 'i', 'j')
  Hier.plot<-Hier.plot[, !names(Hier.plot) %in% drops]
  #reorder the column
  Hier.plot<-Hier.plot[, c(1:7, 9:11, 8, 12:15)]
  return (Hier.plot)
}

#########################################################################################################
#########################################################################################################

#' @rdname Hiermodel
#' @export
Hiergetpi<-function(aedata, hierraw){
  # this function is to get the pit and pic from the output of Hier_history
  # the output is used for Loss function
  # aedata is the output of preprocess
  # hierraw is the output of Hier_history
  # it will output a list with two elements, pit and pic
  # pit contains the columns of SoC, PT for each AE and also the incidence rate
  # for AE in treatment group for each iteraction
  # pic contains the columns of SoC, PT for each AE and also the incidence rate
  # for AE in control group for each iteraction

  # first get \gamma_{bj} and \theta_{bj}
  sim.gamma <- sim.theta <- matrix(0,nrow(aedata),nrow(hierraw))

  for (i in 1:nrow(aedata)){
    ind1 <- aedata$b[i]; ind2 <- aedata$j[i]
    sim.gamma[i,] <- as.numeric(as.character(hierraw[[paste0('gamma.',ind1,'.',ind2,'.')]]))
    sim.theta[i,] <- as.numeric(as.character(hierraw[[paste0('theta.',ind1,'.',ind2,'.')]]))
  }

  # create pit and pic
  # pic=exp(sim.gamma)/(1+exp(sim.gamma))
  # pit=exp(sim.gamma+sim.theta)/(1+exp(sim.gamma+sim.theta))
  pic0<-exp(sim.gamma)/(1+exp(sim.gamma))
  pit0<-exp(sim.theta+sim.gamma)/(1+exp(sim.theta+sim.gamma))

  # combine SoC and PT together with incidence rate
  pit<-cbind(as.character(aedata$AEBODSYS), as.character(aedata$AEDECOD), pit0)
  pic<-cbind(as.character(aedata$AEBODSYS), as.character(aedata$AEDECOD), pic0)
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

#' @rdname Hiermodel
#' @export
Hierplot<-function(hierdata, ptnum=10, param="risk difference", OR_xlim=c(0,5) ){
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
               , legend.background = element_rect(fill="lightblue"))

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
    p<-p+ggtitle(paste0("Top ", ptnum, " AE of median odds ratio \nplotted with 95% credible interval"))
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
               , legend.background = element_rect(fill="lightblue"))

    return(p)
  }
}




#########################################################################################################
#########################################################################################################

#' @rdname Hiermodel
#' @export
Hiertable<-function(hierdata, ptnum=10, param="risk difference" ){
  # this function takes the same parameter as function Hierplot
  # it returns the detailed information of AEs plotted in Hierplot

  library(data.table)
  inputdata<-hierdata

  if (param=="risk difference"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "Diff_mean"], decreasing=TRUE), ]
    test<-head(test, ptnum)
  }

  if (param=="odds ratio"){

    # first to get the top 10 AEs
    test<-inputdata[order(inputdata[, "OR_median"], decreasing=TRUE), ]
    test<-head(test, ptnum)
  }
  return(test)
}



