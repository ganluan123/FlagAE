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
#' The loss is calcuated in following way: first the subjects original AE dataset (output of \code{\link{preprocess2}}) is randomly evenly
#' divided k independent subparts. For each subpart, use this subpart as the testing dataset and use the rest of the whole dataset as the
#' training dataset. Model is trained with the training dataset and then loss is calculated for the testing dataset and training dataset.
#' Repeat this for each subpart and take the average of the testing loss and training loss from each subpart as the final loss. \cr
#'
#' \strong{\code{Lossfun}} takes the AE dataset and fitted incidence as parameters and calculate the loss based on the loss function above.\cr
#'
#' \strong{\code{kfdpar}} first calls function \code{preprocess} to process the data and produce a temporary dataset
#' and also calls function \code{preprocess2} to process the data to get the whole AE dataset.
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
#' @param aedata output from function \code{\link{preprocess2}}
#' @param PI output from function \code{\link{Hiergetpi}} or \code{\link{Isinggetpi}}
#' @param TreatCol a string, used for function \code{\link{preprocess}} and \code{\link{preprocess2}}
#' @param drug a string, used for function \code{\link{preprocess}} and \code{\link{preprocess2}}
#' @inheritParams preprocess2
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
#' AEdata<-preprocess2(adsl=ADSL, adae=ADAE, TreatCol="TREATMENT", drug="xyz")
#' AELIST<-kfdpar(ADSL, ADAE, TreatCol="TREATMENT", drug="xyz", k=5)
#'
#' # Bayesian Hierarchical Model
#' INITS1<-list(mu.gamma.0=0.1, tau.gamma.0=0.1, mu.theta.0=0.1, tau.theta.0=0.1, alpha.pi=2, beta.pi=2)
#' INITS2<-list(mu.gamma.0=1, tau.gamma.0=1, mu.theta.0=1, tau.theta.0=1, alpha.pi=10, beta.pi=10)
#' INITS <- list(INITS1,INITS2)
#' HIERRAW<-Hier_history(aedata=AEdata, inits=INITS, n_burn=1000, n_iter=1000, thin=20, n_adapt=1000, n_chain=2)
#' HIERPI<-Hiergetpi(aedata=AEdata, hierraw=HIERRAW)
#' loss_1<-Lossfun(aedata=AEdata, PI=HIERPI)
#' LOSSHIER<-CVhier(AElist=AELIST, inits=INITS, n_burn=1000, n_iter=1000, thin=20, n_adapt=1000, n_chain=2)
#' LOSSHIER$trainloss # train loss
#' LOSSHIER$testloss # test loss
#'
#' # Bayesian model with Ising prior
#' RHO<-rep(1,dim(AEdata)[1])
#' THETA<-0.02
#' SIM<-c(5000,1000,20)
#' BETA.AB<-c(0.25, 0.75)
#' ISINGRAW<-Ising_history(aedata = AEdata, beta.ab = BETA.AB, rho = RHO, theta = THETA, sim = SIM)
#' ISINGPI<-Isinggetpi(aedata = AEdata, isingraw=ISINGRAW)
#' loss_2<-Lossfun(aedata=AEdata, PI=ISINGPI)
#' LOSSISING<-CVising(AElist=AELIST, beta.ab = BETA.AB, rho = RHO, theta = THETA, sim = SIM)
#' LOSSISING$trainloss # train loss
#' LOSSISING$testloss # test loss
#' }
#'
#' @seealso
#' \code{\link{preprocess}}, \code{\link{preprocess2}}, \code{\link{Hier}}, \code{\link{Ising}}, \code{\link{Isinggetpi}},
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
  # aedata is the output from function preprocess2
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

kfdpar<-function (adsl, adae, TreatCol, drug, k){
  # the first four parameters of this function are the same as in function preprocess
  # it will first call function preprocess to process the data and produce a temporary dataset
  # and also call function preprocess2 to process the data to get the whole AE dataset
  # then this temporary dataset will be  randomly divided into k equal parts
  # this function will generate a list with k elements with each element is a also a list
  # a list contains two elements, named traindf and testdf
  # traindf is used to train the model and testdf is usesd to calcualte the loss
  # this result is going to be used for further crossvalidation to calculate loss


  # get preprocessed dataset
  tempdf<-preprocess(adsl=adsl, adae=adae, TreatCol=TreatCol, drug=drug)
  tempaedata<-preprocess2(adsl=adsl, adae=adae, TreatCol=TreatCol, drug=drug)

  # get the the set of unique usubjid and shuffle it
  IDvect<-unique(tempdf$USUBJID)
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

  # get the AE data with ID inside each element if IDlist
  # and convert into the dataset that can be used for further analysis
  AElist<-list()
  for (i in 1:k){
    # first get the AEs with ID inside the IDlist[[i]]
    Cdat0<-tempdf[tempdf$USUBJID %in% IDlist[[i]], ]

    # the following code is directly copy from function preprocess2

    ################### from preprocess2 ###########################

    # first we need to split Cdat into two dataset
    # one dataset contains only subjects with AE, named cdat
    # the other dataset contains all unique subjects, named adslid
    Cdat<-Cdat0[!is.na(Cdat0$AEBODSYS),]
    # to get adslid dataset, first only take 3 columes
    Adslid<-Cdat0[c("USUBJID", TreatCol, "trt")]
    # then remove the duplicated terms
    library(dplyr)
    Adslid<-Adslid %>% distinct()


    # sum up the number of entries for each combination of AEBODSYS, AEDECOD, and trt
    Tdat <- count(Cdat,AEBODSYS,AEDECOD,trt)
    library(tidyr)
    # get the number of entries for each combination of AEBODSYS and AEDECOD for trt=1 and trt=0, respectively
    Tdat2 <- spread(Tdat, trt, n)

    # for some extreme situation, all the patients are from one group
    # either all from treatment group or all from control group
    if (!( "0" %in% names(Tdat2))) Tdat2["0"]<-NA
    if (!( "1" %in% names(Tdat2))) Tdat2["1"]<-NA

    # change the column names for Tdat2 with trt=1 be AEc, which means # of AEs in control groups
    # and trt=0 be AEt, which means # of AEs in treatment group
    colnames(Tdat2)[colnames(Tdat2) == '0'] <- 'AEc'
    colnames(Tdat2)[colnames(Tdat2) == '1'] <- 'AEt'

    # replave NA by zero
    Tdat2[is.na(Tdat2)] <- 0

    # get the total number of subjects in control and treatment group
    # regardless of having AE or not
    # record the information as Nc and Nt respectively as two columns in dataset
    Tdat2$Nc<-nrow(Adslid[Adslid$trt==0,])
    Tdat2$Nt<-nrow(Adslid[Adslid$trt==1,])

    ###################end of code from preprocess2#########################################

    # A protential problem is that the number of AE in differnt parts may not be the same,
    # and we cannot calculate loss in this case
    # thus we need to make sure all AE shown up in Tdat2, and let AEc and AEt be 0
    # for AE not shown

    # first delete the columns of AEc, AEt, Nc, Nt in tempaedata
    parttempaedata<-tempaedata[, c(1,2,7,8,9)]
    # merge Tdat2 and parttempaedata, and keep all the rows in parttempaedata
    # this is to make sure that all AEs with shown up the final dataset
    Tdat2_tempaedata<-merge(Tdat2, parttempaedata, all.y = TRUE)

    # replace NA in cloumn Nc and Nt with the correct Nc and Nt
    Tdat2_tempaedata$Nc<-Tdat2$Nc[1]
    Tdat2_tempaedata$Nt<-Tdat2$Nt[1]

    # replace all NA in AEc and AEt with 0
    Tdat2_tempaedata[is.na(Tdat2_tempaedata)]<-0

    # Tdat2_tempaedata is the testing dataset
    testset<-Tdat2_tempaedata

    # sort testset by b and j
    testset<-testset[order(testset$b, testset$j),]

    # get training dataset by substract AEc, AEt, Nc, Nt of Tdat2 from tempaedata
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

  return (AElist)
}



#########################################################################################################
#########################################################################################################

#' @rdname crossvalidation
#'
#' @export

CVhier<-function(AElist, inits, n_burn, n_iter, thin, n_adapt, n_chain){
  # parameter AElist is the output from function kfoldpartition
  # the rest parameters are the same as paramters in function Hier_history with the same name
  # this function will calculate the train loss and test loss for each partition of the dataset
  # and return the sum as the final train and test loss

  trainloss<-0
  testloss<-0

  for (i in 1:length(AElist)){
    train<-AElist[[i]]$train
    test<-AElist[[i]]$test

    # train the model
    train_hier<-Hier_history(train, inits, n_burn, n_iter, thin, n_adapt, n_chain)

    # get pi
    train_hiergetpi<-Hiergetpi(aedata=train, hierraw = train_hier)

    # train loss and test loss
    temp_trainloss<-Lossfun(aedata=train, PI=train_hiergetpi)
    temp_testloss<-Lossfun(aedata=test, PI=train_hiergetpi)

    trainloss<-trainloss+temp_trainloss
    testloss<-testloss+temp_testloss
  }
  return(list(trainloss=trainloss/length(AElist), testloss=testloss/length(AElist)))
}

#########################################################################################################
#########################################################################################################

#' @rdname crossvalidation
#'
#' @export

CVising<-function(AElist, beta.ab, rho, theta, sim){
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
    train_ising<-Ising_history(train, beta.ab = beta.ab, rho = rho, theta = theta, sim = sim)

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
