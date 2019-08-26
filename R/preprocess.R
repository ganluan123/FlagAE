#' @name PREPROCESS
#'
#' @title  Raw data process
#'
#' @description
#' \code{preprocess} reorganizes the datasets by adverse events (AE), so that it can be used for further analysis.\cr
#' \code{PREplot} plot the risk difference or odds ratio calculated from Raw data.
#'
#'
#'
#' @details
#' \code{preprocess} sums up the number of patients in control and treatment group.
#' Also this function calculates the number of patients
#' experiencing the AE for each AE in control and treatment group. It also contains columns, "b", "i", "j" for further analysis.\cr
#' \code{PREplot} calculate the risk difference or odds ratio from Raw data and plot top AEs with highest raw risk difference or with highest
#' odds ratio. User can select which summary statistics('risk difference' or 'odds ratio') the plot the based on and also the number
#' of adverse events to plot out.
#'
#' @return
#' \code{preprocess} returns a dataframe with following columnss:\cr
#' \emph{AEBODSYS}: Body System or SoC for AE \cr
#' \emph{AEDECOD}: Peferred Term (PT) for AE \cr
#' \emph{AEc}: number of patients in control group having each specific AE \cr
#' \emph{AEt}: number of patients in treatment group having each specific AE \cr
#' \emph{Nc}: number of total patients in control group \cr
#' \emph{Nt}: number of total patients in treatmetn group \cr
#' \emph{b}: integer represents each Soc \cr
#' \emph{i}: integer represents each PT \cr
#' \emph{j}: order of PT in each SoC \cr
#' \emph{Raw_Risk_Diff}: Risk difference between treatment and control group \cr
#' \emph{Raw_OR}: Odds ratio of the risk treatment/control
#'
#' @param adsl subject level analysis dataset, it is a .csv file, it has to contain at least two columns, "USUBJID" and "TRTCTR", "TRTCTR"
#' is the indicator for treatment and control group. TRTCTR=1 for treatment group and TRTCTR=0 for control group.
#' @param adae adverse event analysis dataset, it is a .csv file, it has to contain at least three columns, "USUBJID", "AEBODSYS",and "AEDECOD"
#' @param aedata the output dataset from function \code{preprocess}
#' @param ptnum an integer, number of PTs to plot, with default be 50
#' @param param a string, summary statistics based for plotting, either 'risk difference' or 'odds ratio', with default be 'risk difference'
#'
#' @examples
#' \dontrun{
#' data(ADAE)
#' data(ADSL)
#' AEdata<-preprocess(adsl=ADSL, adae=ADAE)
#' PREplot(aedata=AEdata)
#' # user can use a very big number(bigger than total PTs in dataset) to plot out all the PTs
#' PREplot(aedata=AEdata, ptnum=50000, param='odds ratio')
#' }
#'
#' @note
#' Make sure that adsl and adae contain the required columns with exact the same column names as listed in "parameters" section.
#'
#' @export
preprocess<-function(adsl, adae){
  ### this function will take 2 parameters
  ### adsl: subject level analysis dataset, it is a .csv file
  ### adae: adverse event analysis dataset, it is a .csv file

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

  # sum up the number of entries for each combination of AEBODSYS, AEDECOD, and TRTCTR
  library(dplyr)
  Tdat <- count(AE,AEBODSYS,AEDECOD,TRTCTR)
  library(tidyr)
  # get the number of entries for each combination of AEBODSYS and AEDECOD for TRTCTR=1 and TRTCTR=0, respectively
  Tdat2 <- spread(Tdat, TRTCTR, n)

  # change the column names for Tdat2 with TRTCTR=1 be AEc, which means # of AEs in control groups
  # and TRTCTR=0 be AEt, which means # of AEs in treatment group
  colnames(Tdat2)[colnames(Tdat2) == '1'] <- 'AEt'
  colnames(Tdat2)[colnames(Tdat2) == '0'] <- 'AEc'

  # replave NA by zero
  Tdat2[is.na(Tdat2)] <- 0

  # get the total number of subjects in control and treatment group
  # regardless of having AE or not
  # record the information as Nc and Nt respectively as two columns in dataset
  Tdat2$Nc<-nrow(SL[SL$TRTCTR==0,])
  Tdat2$Nt<-nrow(SL[SL$TRTCTR==1,])

  # for further analysis get 3 more columns for the dataset
  # b is a column of integers with each integer represent one Soc
  # i is a column of integers with each integer represent one PT
  # j is a column of integers with each integer which is the order of PT in each SoC
  Tdat2['b']<-as.integer(factor(Tdat2$AEBODSYS))
  Tdat2['i']<-as.integer(factor(Tdat2$AEDECOD))
  Tdat2<-Tdat2 %>% group_by(b) %>% mutate(j = as.integer(factor(AEDECOD)))

  # add two columns to dataset,
  # Raw_Risk_Diff, risk difference directly calculated from data
  # Raw_OR, Odds ratio directly calculated from data
  pi_treat<-Tdat2$AEt/Tdat2$Nt
  pi_ctrl<-Tdat2$AEc/Tdat2$Nc
  Tdat2$Raw_Risk_Diff<-pi_treat-pi_ctrl
  odds_treat<-pi_treat/(1-pi_treat)
  odds_ctrl<-pi_ctrl/(1-pi_ctrl)
  Tdat2$Raw_OR<-odds_treat/odds_ctrl

  # make sure Tdat2 is order by b and j for further analysis
  Tdat2<-Tdat2[order(Tdat2$b, Tdat2$j), ]

  # return the dataframe
  return (Tdat2)
}


#########################################################################################################
#########################################################################################################

#' @rdname PREPROCESS
#' @export
PREplot<-function(aedata, ptnum=50, param='risk difference'){
  if(param=='risk difference'){
    library(data.table)
    library(ggplot2)
    # reoder aedata by model based risk difference
    aedata<-aedata[order(-aedata$Raw_Risk_Diff), ]
    if(ptnum>(dim(aedata)[1])) ptnum<-dim(aedata)[1]
    test<-copy(aedata[1:ptnum, ])
    test$x<-1:ptnum
    p<-ggplot(test, aes(x))+geom_point(aes(y=Raw_Risk_Diff),colour='#3591d1')
    p<-p+theme(legend.title=element_blank())
    p<-p+ggtitle(paste0('Risk difference from raw data'))
    p<-p+xlab("PT")+ylab('Risk difference')

  }

  if(param=='odds ratio'){
    library(data.table)
    library(ggplot2)
    # reoder aedata by model based odds ratio
    aedata<-aedata[order(-aedata$Raw_OR), ]
    if(ptnum>(dim(aedata)[1])) ptnum<-dim(aedata)[1]
    test<-copy(aedata[1:ptnum, ])
    test$x<-1:ptnum
    p<-ggplot(test, aes(x))+geom_point(aes(y=Raw_OR),colour='#3591d1')
    p<-p+theme(legend.title=element_blank())
    p<-p+ggtitle(paste0('Odds ratio from raw data'))
    p<-p+xlab("PT")+ylab('Odds ratio')
  }
  p<-p + theme(axis.title.x = element_text(size=13)) + theme(axis.title.y = element_text(size=13))
  p<-p + theme(plot.title = element_text(size=15, hjust=0.5))
  p<-p + theme(axis.text.x=element_text(size=15)) + theme(axis.text.y=element_text(size=15))
  return(p)
}
