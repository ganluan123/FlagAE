Hier_history<- function(aedata, n_burn, n_iter, thin, n_adapt, n_chain, alpha.gamma=3, beta.gamma=1,
                        alpha.theta=3, beta.theta=1, mu.gamma.0.0=0, tau.gamma.0.0=10, alpha.gamma.0.0=3,
                        beta.gamma.0.0=1, lambda.alpha=0.1, lambda.beta=0.1, mu.theta.0.0=0, tau.theta.0.0=10,
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
  mu.gamma.0 ~ dnorm(mu.gamma.0.0, 1/tau.gamma.0.0)
  tau.gamma.0 ~ dgamma(alpha.gamma.0.0, beta.gamma.0.0)

  # hyperpriors for thetas;
  mu.theta.0 ~ dnorm(mu.theta.0.0, 1/tau.theta.0.0)
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
    temp.fit <- jags.model(textConnection(model.binom),data=data,n.chains=1,
                           inits=list(.RNG.name='base::Wichmann-Hill', .RNG.seed=1), n.adapt=n_adapt,quiet=TRUE)
    update(temp.fit, n.iter=n_burn)

    # summary of posterior samples
    temp.param.samples <- coda.samples(temp.fit,param,n.iter=n_iter,thin=thin)
    temp.param.est <- data.frame(as.matrix(temp.param.samples))
    Final.est<-temp.param.est
  }

  return(Final.est)
}



Hier_history2<- function(aedata, n_burn, n_iter, thin, n_adapt, n_chain, alpha.gamma=3, beta.gamma=1,
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









#INITS1<-list(mu.gamma.0=0.1, tau.gamma.0=0.1, mu.theta.0=0.1, tau.theta.0=0.1, alpha.pi=2, beta.pi=2)
#INITS2<-list(mu.gamma.0=1, tau.gamma.0=1, mu.theta.0=1, tau.theta.0=1, alpha.pi=10, beta.pi=10)
#INITS <- list(INITS1,INITS2)
HIERRAW<-Hier_history(aedata=AEdata,  n_burn=1000, n_iter=100, thin=2, n_adapt=1000, n_chain=1)

HIERRAW2<-Hier_history2(aedata=AEdata,  n_burn=1000, n_iter=100, thin=2, n_adapt=1000, n_chain=1)

HIERRAW[,3]-HIERRAW2[,3]
