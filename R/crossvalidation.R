#########################################################################################################
#########################################################################################################

#' @name crossvalidation
#'
#' @title Loss Calculation by Cross Validation
#'
#' @description Function here are to calculate the loss by cross validation for Bayesian hierarchical model (see also \code{\link{Hier}})
#' and Bayesian model with Ising prior (see also \code{\link{Ising}}). This can be used to select the best hyperparameters and to compare
#' two models.
#'
#' @details
#' The loss is calcuated by:
#' \deqn{\sqrt{\sum_{bj} [(Y_{bj}-N_t*t_{bj})^2]}/N_t + \sqrt{\sum_{bj} [(X_{bj}-N_c*c_{bj})^2]}/N_c}
#' Here b=1,..., B and j=1, ... , k_b, Y_{bj} and X_{bj} are the number of subjects with
#' an AE with PT j under SOC b in treatment and control groups.
#' N_t and N_c are the number of subjects in treatment and control groups, respectively.
#' t_{bj} and c_{bj} are the model fitted incidence of an AE with PT j under SOC b in treatment and control groups.
#' This formular gives the loss for one interaction/sample, the final loss is the average of loss from all of the interactions/samples.\cr
#'
#' The loss is calcuated in following way: first the subjects original AE dataset (output of \code{\link{preprocess}}) is randomly evenly
#' divided k independent subparts. For each subpart, use this subpart as the testing dataset and use the rest of the whole dataset as the
#' training dataset. Model is trained with the training dataset and then loss is calculated for the testing dataset and training dataset.
#' Repeat this for each subpart and take the average of the testing loss and training loss from each subpart as the final loss. \cr
#'
#' \strong{\code{Lossfun}} takes the AE dataset and fitted incidence as parameters and calculate the loss based on the loss function above.\cr
#'
#' \strong{\code{kfdpar}} first calls function \code{preprocess} to process the data and produce a temporary dataset
#' and also calls function \code{preprocess} to process the data to get the whole AE dataset.
#' Then this temporary dataset will be  randomly divided into k equal subparts. For each subpart,
#' use this subpart as the testing dataset and use the rest of the whole dataset as the
#' training dataset.This function will generate a list with k elements with each element is a also a list
#' a list contains two elements, named traindf and testdf.
#' "traindf" is used to train the model and testdf is usesd to calcualte the loss.
#' The output is going to be used for further crossvalidation to calculate loss. \cr
#'
#' \strong{\code{CVhier}} calculates the loss for Bayesian Hierarchical model.\cr
#'
#' \strong{\code{CVising}} calculates the los for Bayesian model with Ising prior. \cr
#'
#' @inheritParams preprocess
#' @param aedata output from function \code{\link{preprocess}}
#' @param PI output from function \code{\link{Hiergetpi}} or \code{\link{Isinggetpi}}
#' @param k interger, the number of folds used to split the dataset for cross validation
#' @inheritParams Ising
#' @inheritParams Hier
#'
#' @return
#' \strong{\code{Lossfun}} returns the loss for dataset \code{aedata} based on the fitted incidence \code{PI}.\cr
#' \strong{\code{kfdpar}} returns a list with k elements with each element is a also a list,
#' that contains two elements, named traindf and testdf.\cr
#' \strong{\code{CVhier}} returns the final training and testing loss for Bayesian hierarchical model. \cr
#' \strong{\code{CVIsing}} returns the final training and testing loss for Bayesian model with Ising model. \cr
#'
#' @examples
#' \dontrun{
#' data(ADAE)
#' data(ADSL)
#' AEdata<-preprocess(adsl=ADSL, adae=ADAE)
#' AELIST<-kfdpar(ADSL, ADAE, k=5)
#'
#' # Bayesian Hierarchical Model
#' HIERRAW<-Hier_history(aedata=AEdata, n_burn=1000, n_iter=1000, thin=20, n_adapt=1000, n_chain=2)
#' HIERPI<-Hiergetpi(aedata=AEdata, hierraw=HIERRAW)
#' loss_1<-Lossfun(aedata=AEdata, PI=HIERPI)
#' LOSSHIER<-CVhier(AElist=AELIST, n_burn=1000, n_iter=1000, thin=20, n_adapt=1000, n_chain=2)
#' LOSSHIER$trainloss # train loss
#' LOSSHIER$testloss # test loss
#'
#' # Bayesian model with Ising prior
#' RHO<-rep(1,dim(AEdata)[1])
#' ISINGRAW<-Ising_history(aedata = AEdata, n_burn=1000, n_iter=5000, thin=20, alpha=0.5, beta=0.75, alpha.t=0.5, beta.t=0.75,
#'                                    alpha.c=0.25, beta.c=0.75, rho=RHO, theta=0.02)
#' ISINGPI<-Isinggetpi(aedata = AEdata, isingraw=ISINGRAW)
#' loss_2<-Lossfun(aedata=AEdata, PI=ISINGPI)
#'
#' LOSSISING<-CVising(AElist=AELIST, n_burn=100, n_iter=500, thin=20, alpha=0.5, beta=0.75, alpha.t=0.5, beta.t=0.75,
#'                              alpha.c=0.25, beta.c=0.75, rho=RHO, theta=0.02)
#' LOSSISING$trainloss # train loss
#' LOSSISING$testloss # test loss
#' }
#'
#' @seealso
#' \code{\link{preprocess}}, \code{\link{Hier}}, \code{\link{Ising}}, \code{\link{Isinggetpi}},
#' \code{\link{Hiergetpi}}
#'
#' @export

Lossfun<-function (aedata, PI){
  # aedata is the output from function preprocess
  # PI is the output from either function Hiergetpi or Isiinggetpi
  # PI is a list containg two items, pit and pic
  # this function will calculate the loss for dataset aedata,
  # based on the predicated incidence rate for treatment group pit
  # and the predicated incidence rate for control group pic
  # aedata is the output from function preprocess
  # pit dataset must contain the following items (and only contain these items)
  # SoC, PT, incidience rate for treatment group from each iteraction (SoC and PT must be the first two columns)
  # pic dataset must contain the following items
  # SoC, PT, incidience rate for control group from each iteraction (SoC and PT must be the first two columns)
  # the loss is defined as:  $\sqrt{\sum_{bj} [(Y_{bj}-N_t*t_{bj})^2]}/N_t + \sqrt{\sum_{bj} [(X_{bj}-N_c*c_{bj})^2]}/N_c$

  # get pit and pic from PI
  pit<-PI$pit
  pic<-PI$pic
  # sort the dataset so that AE in these dataset are in the same order
  aedata<-aedata[order(aedata$AEBODSYS, aedata$AEDECOD), ]
  pit<-pit[order(pit$SoC, pit$PT),]
  pic<-pic[order(pic$SoC, pic$PT),]

  # it is possible that the number in pit/pic is stored as factor
  # convert it to numeric
  for (i in 3:(dim(pit)[2])){
    pit[,i]<-as.numeric(as.character(pit[,i]))
    pic[,i]<-as.numeric(as.character(pic[,i]))
  }

  # remove SOC and PT, left only predicted incidence rate
  pit2<-pit[, 3:(dim(pit)[2])]
  pic2<-pic[, 3:(dim(pic)[2])]

  # calculate the loss contributed by treatment group
  loss_t<-sqrt(mean(colSums((aedata$AEt - aedata$Nt*pit2)^2)))/aedata$Nt[1]

  # calculate the loss contributed by control group
  loss_c<-sqrt(mean(colSums((aedata$AEc - aedata$Nc*pic2)^2)))/aedata$Nc[1]

  # get the total loss
  loss<-loss_t+loss_c

  return (loss)
}


#########################################################################################################
#########################################################################################################

#' @rdname crossvalidation
#'
#' @export

kfdpar<-function (adsl, adae, k){
  # the first two parameters of this function are the same as in function preprocess
  # k is integer which determine how many folds will the dataset be divided
  # this function will randomly divide the subjects into k equal parts
  # for each part i, this subjects in this part is assigned to testing dataset
  # the rest (k-1) parts are assigned to training dataset
  # this function will generate a list with k elements with each element is a also a list
  # a list contains two elements, named traindf and testdf
  # traindf is used to train the model and testdf is usesd to calcualte the loss
  # this result is going to be used for further crossvalidation to calculate loss

  # get the AEdata for whole dataset
  tempaedata<-preprocess(adsl=adsl, adae=adae)

  ### only extract columns "USUBJID" and "TRTCTR" from adsl
  ### USUBJID: Unique Subject Identifier
  ### TRTCTR: indicator for treatment (TRTCTR=1) or control (TRTCTR=0) group
  SL<-adsl[!is.na(adsl$TRTCTR), ]
  SL<-subset(adsl, select=c("USUBJID","TRTCTR"))
  SL<-unique(SL)


  ### extract columns "USUBJID", "AEBODSYS", and "AEDECOD"
  ### USUBJID: Unique Subject Identifier
  ### AEBODSYS: Body System or Organ Class (SoC)
  ### AEDECOD: Dictionary-Derived Term (Perferred Term, PT)
  aevars<-c("USUBJID","AEBODSYS","AEDECOD")
  AE<-subset(adae, select=aevars)
  AE<-unique(AE)
  ### order AE by AEBODSYS and AEDECOD
  AE<-AE[with(AE, order(AEBODSYS, AEDECOD)), ]


  ### merge SL and AE by usubjid
  ### we will discard the entries in AE that with usubjid not in SL
  IS<-intersect(unique(AE$USUBJID), SL$USUBJID)
  AE<-AE[AE$USUBJID %in% IS, ]

  ### add column "TRTCTR" for AE
  AE$TRTCTR<-0
  for (i in 1:dim(AE)[1]){
    ID<-as.character(AE[i, "USUBJID"])
    trtcode<-SL[SL$USUBJID %in% ID, "TRTCTR"]
    AE[i, "TRTCTR"]<-trtcode
  }

  # get the the set of unique usubjid and shuffle it
  # IDvect<-unique(tempdf$USUBJID)
  IDvect<-unique(SL$USUBJID)
  IDvect<-sample(IDvect)

  # divide IDdf into k parts, if the number of ID is not exactly divisible, put the residual into the last part
  IDnum<-length(IDvect)
  eachnum<-floor(IDnum/k)
  IDlist<-list()
  for (i in 1:(k-1)){
    Index<-(eachnum*(i-1)+1):(eachnum*i)
    IDlist[[i]]<-IDvect[Index]
  }
  IDlist[[k]]<-IDvect[(eachnum*(k-1)+1):IDnum]

  # get the AE data with ID inside each element of IDlist
  # and convert into the dataset that can be used for further analysis
  AElist<-list()
  for (i in 1:k){

    # first get the tempAE and tempSL with USUBJID in IDlist[[i]]
    tempAE<-AE[AE$USUBJID %in% IDlist[[i]], ]
    tempSL<-SL[SL$USUBJID %in% IDlist[[i]], ]

    # sum up the number of entries for each combination of AEBODSYS, AEDECOD, and TRTCTR
    library(dplyr)
    tempTdat <- count(tempAE,AEBODSYS,AEDECOD,TRTCTR)
    library(tidyr)
    # get the number of entries for each combination of AEBODSYS and AEDECOD for TRTCTR=1 and TRTCTR=0, respectively
    tempTdat2 <- spread(tempTdat, TRTCTR, n)

    # change the column names for tempTdat2 with TRTCTR=1 be AEc, which means # of AEs in control groups
    # and TRTCTR=0 be AEt, which means # of AEs in treatment group
    colnames(tempTdat2)[colnames(tempTdat2) == '1'] <- 'AEt'
    colnames(tempTdat2)[colnames(tempTdat2) == '0'] <- 'AEc'

    # replave NA by zero
    tempTdat2[is.na(tempTdat2)] <- 0

    # get the total number of subjects in control and treatment group
    # regardless of having AE or not
    # record the information as Nc and Nt respectively as two columns in dataset
    tempTdat2$Nc<-nrow(tempSL[tempSL$TRTCTR==0,])
    tempTdat2$Nt<-nrow(tempSL[tempSL$TRTCTR==1,])

    # A protential problem is that the number of AE in differnt parts may not be the same,
    # and we cannot calculate loss in this case
    # thus we need to make sure all AE shown up in Tdat2, and let AEc and AEt be 0
    # for AE not shown

    # first delete the columns of AEc, AEt, Nc, Nt in tempaedata
    parttempaedata<-tempaedata[, c("AEBODSYS", "AEDECOD", "b", "i", "j")]
    # merge tempTdat2 and parttempaedata, and keep all the rows in parttempaedata
    # this is to make sure that all AEs with shown up the final dataset
    Tdat2_tempaedata<-merge(tempTdat2, parttempaedata, all.y = TRUE)

    # replace NA in cloumn Nc and Nt with the correct Nc and Nt
    Tdat2_tempaedata$Nc<-tempTdat2$Nc[1]
    Tdat2_tempaedata$Nt<-tempTdat2$Nt[1]

    # replace all NA in AEc and AEt with 0
    Tdat2_tempaedata[is.na(Tdat2_tempaedata)]<-0

    # Tdat2_tempaedata is the testing dataset
    testset<-Tdat2_tempaedata

    # sort testset by b and j
    testset<-testset[order(testset$b, testset$j),]

    # get training dataset by substract AEc, AEt, Nc, Nt of testset from tempaedata
    train<-tempaedata
    #train[,c("AEc", "AEt", "Nc", "Nt")]<-tempaedata[, c("AEc", "AEt", "Nc", "Nt")]-testset[, c("AEc", "AEt", "Nc", "Nt")]
    train[,"AEc"]<-tempaedata[, "AEc"]-testset[, "AEc"]
    train[,"AEt"]<-tempaedata[, "AEt"]-testset[, "AEt"]
    train[,"Nc"]<-tempaedata[, "Nc"]-testset[, "Nc"]
    train[,"Nt"]<-tempaedata[, "Nt"]-testset[, "Nt"]

    # put train and testset toghether in a list
    train_test<-list(train=train, test=testset)

    # add train_test into AElist
    AElist[[i]]<-train_test
  }
  return(AElist)
}

#########################################################################################################
#########################################################################################################

#' @rdname crossvalidation
#'
#' @export

CVhier<-function(AElist, n_burn, n_iter, thin, n_adapt, n_chain,
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

  library(foreach)
  library(doParallel)
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
#########################################################################################################
#########################################################################################################

#' @rdname crossvalidation
#'
#' @export

CVising<-function(AElist, n_burn, n_iter, thin, alpha_=0.25, beta_=0.75, alpha.t=0.25, beta.t=0.75,
                  alpha.c=0.25, beta.c=0.75, rho, theta){
  # parameter AElist is the output from function kfoldpartition
  # the rest parameters are the same as paramters in function Ising_history with the same name
  # this function will calculate the train loss and test loss for each partition of the dataset
  # and return the sum as the final train and test loss
  trainloss<-0
  testloss<-0

  for (i in 1:length(AElist)){
    train<-AElist[[i]]$train
    test<-AElist[[i]]$test

    # train the model
    train_ising<-Ising_history(aedata =train, n_burn=n_burn, n_iter=n_iter, thin=thin,
                               alpha_=alpha_, beta_=beta_, alpha.t=alpha.t, beta.t=beta.t,
                               alpha.c=alpha.c, beta.c=beta.c, rho=rho, theta=theta)

    # get pi
    train_isinggetpi<-Isinggetpi(aedata=train, isingraw = train_ising)

    # train loss and test loss
    temp_trainloss<-Lossfun(aedata=train, PI=train_isinggetpi)
    temp_testloss<-Lossfun(aedata=test, PI=train_isinggetpi)

    trainloss<-trainloss+temp_trainloss
    testloss<-testloss+temp_testloss
  }
  return(list(trainloss=trainloss/length(AElist), testloss=testloss/length(AElist)))
}
