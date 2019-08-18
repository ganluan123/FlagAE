library(FlagAE)

CVhier2<-function(AElist, n_burn, n_iter, thin, n_adapt, n_chain,
                 alpha.gamma=3, beta.gamma=1,
                 alpha.theta=3, beta.theta=1, mu.gamma.0.0=0, tau.gamma.0.0=0.1, alpha.gamma.0.0=3,
                 beta.gamma.0.0=1, lambda.alpha=0.1, lambda.beta=0.1, mu.theta.0.0=0, tau.theta.0.0=0.1,
                 alpha.theta.0.0=3, beta.theta.0.0=1){
  # parameter AElist is the output from function kfoldpartition
  # the rest parameters are the same as paramters in function Hier_history with the same name
  # this function will calculate the train loss and test loss for each partition of the dataset
  # and return the sum as the final train and test loss

  trainloss<-0
  testloss<-0


  cores=detectCores()
  cl<-makeCluster(cores[1]-1)
  registerDoParallel(cl)


  LOSSlist<-foreach (i = 1:length(AElist), .combine=rbind) %dopar% {
    train<-AElist[[i]]$train
    test<-AElist[[i]]$test
    library(FlagAE)

    # train the model
    train_hier<-Hier_history(aedata=train, n_burn=n_burn, n_iter=n_iter, thin=thin, n_adapt=n_adapt, n_chain=n_chain, alpha.gamma=alpha.gamma, beta.gamma=beta.gamma,
                             alpha.theta=alpha.theta, beta.theta=beta.theta, mu.gamma.0.0=mu.gamma.0.0, tau.gamma.0.0=tau.gamma.0.0,
                             alpha.gamma.0.0=alpha.gamma.0.0, beta.gamma.0.0=beta.gamma.0.0, lambda.alpha=lambda.alpha,
                             lambda.beta=lambda.beta, mu.theta.0.0=mu.theta.0.0, tau.theta.0.0=tau.theta.0.0,
                             alpha.theta.0.0=alpha.theta.0.0, beta.theta.0.0=beta.theta.0.0)

    # get pi
    train_hiergetpi<-Hiergetpi(aedata=train, hierraw = train_hier)

    # train loss and test loss
    temp_trainloss<-Lossfun(aedata=train, PI=train_hiergetpi)
    temp_testloss<-Lossfun(aedata=test, PI=train_hiergetpi)

    c(temp_trainloss, temp_testloss)
  }
  stopCluster(cl)

  return(list(trainloss=mean(LOSSlist[,1]), testloss=mean(LOSSlist[,2])))
}


### compare running time
Start<-Sys.time()
LOSSHIER<-CVhier(AElist=AELIST, n_burn=2000, n_iter=5000, thin=20, n_adapt=1000, n_chain=2)
Sys.time()-Start


Start<-Sys.time()
LOSSHIER2<-CVhier2(AElist=AELIST, n_burn=2000, n_iter=5000, thin=20, n_adapt=1000, n_chain=2)
Sys.time()-Start

