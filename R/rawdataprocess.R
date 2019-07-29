#' @name rawdataprocess
#'
#' @title  Raw data process
#'
#' @description
#' \code{preprocess} divides the subjects(patients) in adsl into two groups
#' by creating a new column 'trt' with trt=1 for treatment group and trt=0 for control group.
#' Also this function only takes the columns related to further analysis. \cr
#' \code{preprocess2} first execute function \code{preprocess}
#' and reorganize the dataset by adverse events (AE), so that it can be used for further analysis.\cr
#'
#' \code{preprocess3} is same as \code{preprocess} except that \code{preprocess3} takes 'treatcode' as
#' input instead of 'drug'. Like 'drug' in \code{preprocess}, 'treatcode' here is also used to determine
#' whether one subject is belong to treatment or control group. If the value under column TreatCol is
#' an element of 'treatcode' then this subject belong to treatment group (trt=1) otherwise it belongs
#' control group (trt=0). \cr
#'
#' \code{preprocess4} is same as \code{preprocess2} except it first exectue function \code{preprocess3}. \cr
#'
#' \code{preprocess3} and \code{preprocess4} provide more flexiblity for processing the raw data.
#'
#'
#' @details
#' Parameter "TreatCol" is the column based on which the division is performed,
#' paramter "drug" is to determine the value of trt,
#' if the value under TreatCol contains drug then trt=1, otherwise trt=0.\cr
#' \code{preprocess2} sums up the number of patients in control and treatment group.
#' Also this function calculates the number of patients in control and treatment group
#' experiencing the AE for each AE. It also contains columns, "b", "i", "j" for further analysis.
#'
#' @return
#' \code{preprocess} returns a dataframe with following columns:\cr
#' \emph{USUBJID}: the unique subject id for all subjects in the trail with or without AE \cr
#' \emph{TreatCol}: the column based on which the treatment/control division is performed \cr
#' \emph{TreatCol}: the column based on which the treatment/control division is performed\cr
#' \emph{AEBODSYS}: Body System or SoC for AE\cr
#' \emph{AEDECOD}: Peferred Term (PT) for AE (if a subject did not have any AE, AEBODSYS and AEDECOD are NA)\cr
#' \emph{trt}: indication of treatment group (1) or control group (0)\cr
#'
#' \code{preprocess2} returns a dataframe with following columnss:\cr
#' \emph{AEBODSYS}: Body System or SoC for AE \cr
#' \emph{AEDECOD}: Peferred Term (PT) for AE \cr
#' \emph{AEc}: number of patients in control group having each specific AE \cr
#' \emph{AEt}: number of patients in treatment group having each specific AE \cr
#' \emph{Nc}: number of total patients in control group \cr
#' \emph{Nt}: number of total patients in treatmetn group \cr
#' \emph{b}: integer represents each Soc \cr
#' \emph{i}: integer represents each PT \cr
#' \emph{j}: order of PT in each SoC \cr
#'
#' @param adsl subject level analysis dataset, it is a .sas7bdat file
#' @param adae adverse event analysis dataset, it is a .sas7bdat file
#' @param TreatCol a string
#' @param drug a string
#' @param treatcode a vector of strings
#'
#' @examples
#' \dontrun{
#' data(ADAE)
#' data(ADSL)
#' preprocess(adsl=ADSL, adae=ADAE, TreatCol="TREATMENT", drug="xyz")
#' preprocess2(adsl=ADSL, adae=ADAE, TreatCol="TREATMENT", drug="xyz")
#' preprocess3(adsl=ADSL, adae=ADAE, TreatCol="TREATMENT", treatcode=c("xyz", "xyz2"))
#' preprocess4(adsl=ADSL, adae=ADAE, TreatCol="TREATMENT", treatcode=c("xyz", "xyz2"))
#'
#' }
#'
#' @note
#' Both functions only select the patients that at least took one dose of the treatment
#' and only select the adverse events occurs within 30 days after taking the first dose. \cr
#' \code{preprocess} summarizes the information about AE, and its output can be used
#' for exploratory analysis. \cr
#' The outout of \code{preprocess2} is organized by column b and column j for further analysis.
#' @export
preprocess<-function(adsl, adae, TreatCol, drug){
  ### this function will take 4 parameters, adsl, adae, TreatCol, drug
  ### adsl: subject level analysis dataset, it is a .sas7bdat file
  ### adae: adverse event analysis dataset, it is a .sas7bdat file
  ### both TreatCol and drug are strings
  ### this function will divide the subject(patient) in adsl into two groups
  ### by creating a new column 'trt' with trt=1 for treatment group
  ### and trt=0 for control group
  ### TreatCol is the column based on which the division is performed
  ### drug is to determine the value of trt, if the value under TreatCol
  ### contains drug then trt=1, otherwise trt=0
  ### this function will return a dataset containg the following columns:
  ### USUBJID: the unique subject id for all subjects in the trail with or without AE
  ### TreatCol: the column based on which the treatment/control division is performed
  ### AEBODSYS: Body System or SoC for AE
  ### AEDECOD: Peferred Term (PT) for AE
  ### if a subject did not have any AE, AEBODSYS and AEDECOD are NA
  ### trt: indication of treatment group (1) or control group (0)


  ### select TEAE
  ### select the adverse events occurs within 30 days after taking the first dose
  adae<-adae[adae$TRTEMFL=='Y',]


  ### select treated subjects
  ### select the patients that at least took one dose of the treatment
  adsl<-adsl[adsl$SAFFL == 'Y',]


  ### only extract columns "USUBJID" and TreatCol from adsl
  ### USUBJID: Unique Subject Identifier
  ### TRT01A: Actual Treatment for Period 01
  adsl1<-subset(adsl, select=c("USUBJID",TreatCol))


  ### extract columns "USUBJID", "AEBODSYS", and "AEDECOD"
  ### USUBJID: Unique Subject Identifier
  ### AEBODSYS: Body System or Organ Class (SoC)
  ### AEDECOD: Dictionary-Derived Term (Perferred Term, PT)
  aevars<-c("USUBJID","AEBODSYS","AEDECOD")
  ae1<-subset(adae, select=aevars)
  ae1<-unique(ae1)
  ### order ae1 by AEBODSYS and AEDECOD
  ae1<-ae1[with(ae1, order(AEBODSYS, AEDECOD)), ]


  ### create a new varible trt to indicate
  ### for the indication of treatment group or control group
  ### if drug is contained in TreatCol for one subject, then
  ### trt=1, otherwise trt=0
  #Treatment<-adsl1[TreatCol]
  Treatment<-subset(adsl, select=TreatCol)
  trt<-rep(0, dim(adsl1)[1])
  for (i in 1:dim(adsl1)[1]){
    if (grepl(drug,Treatment[i,1])) trt[i]<-1
  }
  adsl1$trt<-trt


  ### merge adsl1 and ae1 by usubjid
  ### we will discard the entries in ae1 that with usubjid not in adsl1
  IS<-intersect(unique(ae1$USUBJID), adsl1$USUBJID)
  ae1<-ae1[ae1$USUBJID %in% IS, ]

  ### merge ae1 and adsl1
  ### with subject that has no adverse events have NA for AEBODSYS and AEDECOD
  cdat <- merge(adsl1,ae1,by="USUBJID", all=T)
  cdat <- unique(cdat)

  ### return cdat
  return (cdat)
}


#########################################################################################################
#########################################################################################################

#' @rdname rawdataprocess
#'
#' @export
preprocess2<-function (adsl, adae, TreatCol, drug){
  # this function takes the same parameters as preprocess
  # and this function do further data preparation based on the result from function preprocess
  # this function will return a dataframe with following columns
  # AEBODSYS: Body System or SoC for AE
  # AEDECOD: Peferred Term (PT) for AE
  # AEc: number of patients in control group having each specific AE
  # AEt: number of patients in treatment group having each specific AE
  # Nc: number of total patients in control group
  # Nt: number of total patients in treatmetn group
  # b: integer represents each Soc
  # i: integer represents each PT
  # j: order of PT in each SoC

  Cdat0<-preprocess(adsl=adsl, adae=adae, TreatCol=TreatCol, drug=drug)
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

  # change the column names for Tdat2 with trt=1 be AEc, which means # of AEs in control groups
  # and trt=0 be AEt, which means # of AEs in treatment group
  colnames(Tdat2)[colnames(Tdat2) == '1'] <- 'AEt'
  colnames(Tdat2)[colnames(Tdat2) == '0'] <- 'AEc'

  # replave NA by zero
  Tdat2[is.na(Tdat2)] <- 0

  # get the total number of subjects in control and treatment group
  # regardless of having AE or not
  # record the information as Nc and Nt respectively as two columns in dataset
  Tdat2$Nc<-nrow(Adslid[Adslid$trt==0,])
  Tdat2$Nt<-nrow(Adslid[Adslid$trt==1,])

  # for further analysis get 3 more columns for the dataset
  # b is a column of integers with each integer represent one Soc
  # i is a column of integers with each integer represent one PT
  # j is a column of integers with each integer which is the order of PT in each SoC
  Tdat2['b']<-as.integer(factor(Tdat2$AEBODSYS))
  Tdat2['i']<-as.integer(factor(Tdat2$AEDECOD))
  Tdat2<-Tdat2 %>% group_by(b) %>% mutate(j = as.integer(factor(AEDECOD)))

  # make sure Tdat2 is order by b and j for further analysis
  Tdat2<-Tdat2[order(Tdat2$b, Tdat2$j), ]

  # return the dataframe
  return (Tdat2)
}


#########################################################################################################
#########################################################################################################


#' @rdname rawdataprocess
#'
#' @export

preprocess3<-function(adsl, adae, TreatCol, treatcode){
  ### this function will take 4 parameters, adsl, adae, TreatCol, treatcode
  ### adsl: subject level analysis dataset, it is a .sas7bdat file
  ### adae: adverse event analysis dataset, it is a .sas7bdat file
  ### TreatCol: a string
  ### treatcode: a vector of strings
  ### this function will divide the subject(patient) in adsl into two groups
  ### by creating a new column 'trt' with trt=1 for treatment group
  ### and trt=0 for control group
  ### TreatCol is the column based on which the division is performed
  ### treatcode, a vectors of strings,  is to determine the value of trt,
  ### if the value under TreatCol is inside treatcode, then trt=1, otherwise trt=0
  ### this function will return a dataset containg the following columns:
  ### USUBJID: the unique subject id for all subjects in the trail with or without AE
  ### TreatCol: the column based on which the treatment/control division is performed
  ### AEBODSYS: Body System or SoC for AE
  ### AEDECOD: Peferred Term (PT) for AE
  ### if a subject did not have any AE, AEBODSYS and AEDECOD are NA
  ### trt: indication of treatment group (1) or control group (0)


  ### select TEAE
  ### select the adverse events occurs within 30 days after taking the first dose
  adae<-adae[adae$TRTEMFL=='Y',]


  ### select treated subjects
  ### select the patients that at least took one dose of the treatment
  adsl<-adsl[adsl$SAFFL == 'Y',]


  ### only extract columns "USUBJID" and TreatCol from adsl
  ### USUBJID: Unique Subject Identifier
  ### TRT01A: Actual Treatment for Period 01
  adsl1<-subset(adsl, select=c("USUBJID",TreatCol))


  ### extract columns "USUBJID", "AEBODSYS", and "AEDECOD"
  ### USUBJID: Unique Subject Identifier
  ### AEBODSYS: Body System or Organ Class (SoC)
  ### AEDECOD: Dictionary-Derived Term (Perferred Term, PT)
  aevars<-c("USUBJID","AEBODSYS","AEDECOD")
  ae1<-subset(adae, select=aevars)
  ae1<-unique(ae1)
  ### order ae1 by AEBODSYS and AEDECOD
  ae1<-ae1[with(ae1, order(AEBODSYS, AEDECOD)), ]


  ### create a new varible trt to indicate
  ### for the indication of treatment group or control group
  ### if TreatCol for one subject is contained in treatcode, then
  ### trt=1, otherwise trt=0
  #Treatment<-adsl1[TreatCol]
  Treatment<-subset(adsl, select=TreatCol)
  trt<-rep(0, dim(adsl1)[1])
  for (i in 1:dim(adsl1)[1]){
    if (Treatment[i,1] %in% treatcode) trt[i]<-1
  }
  adsl1$trt<-trt


  ### merge adsl1 and ae1 by usubjid
  ### we will discard the entries in ae1 that with usubjid not in adsl1
  IS<-intersect(unique(ae1$USUBJID), adsl1$USUBJID)
  ae1<-ae1[ae1$USUBJID %in% IS, ]

  ### merge ae1 and adsl1
  ### with subject that has no adverse events have NA for AEBODSYS and AEDECOD
  cdat <- merge(adsl1,ae1,by="USUBJID", all=T)
  cdat <- unique(cdat)

  ### return cdat
  return (cdat)
}


#########################################################################################################
#########################################################################################################

#' @rdname rawdataprocess
#'
#' @export
preprocess4<-function (adsl, adae, TreatCol, treatcode){
  # this function takes the same parameters as preprocess3
  # and this function do further data preparation based on the result from function preprocess
  # this function will return a dataframe with following columns
  # AEBODSYS: Body System or SoC for AE
  # AEDECOD: Peferred Term (PT) for AE
  # AEc: number of patients in control group having each specific AE
  # AEt: number of patients in treatment group having each specific AE
  # Nc: number of total patients in control group
  # Nt: number of total patients in treatmetn group
  # b: integer represents each Soc
  # i: integer represents each PT
  # j: order of PT in each SoC

  Cdat0<-preprocess3(adsl=adsl, adae=adae, TreatCol=TreatCol, treatcode=treatcode)
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

  # change the column names for Tdat2 with trt=1 be AEc, which means # of AEs in control groups
  # and trt=0 be AEt, which means # of AEs in treatment group
  colnames(Tdat2)[colnames(Tdat2) == '1'] <- 'AEt'
  colnames(Tdat2)[colnames(Tdat2) == '0'] <- 'AEc'

  # replave NA by zero
  Tdat2[is.na(Tdat2)] <- 0

  # get the total number of subjects in control and treatment group
  # regardless of having AE or not
  # record the information as Nc and Nt respectively as two columns in dataset
  Tdat2$Nc<-nrow(Adslid[Adslid$trt==0,])
  Tdat2$Nt<-nrow(Adslid[Adslid$trt==1,])

  # for further analysis get 3 more columns for the dataset
  # b is a column of integers with each integer represent one Soc
  # i is a column of integers with each integer represent one PT
  # j is a column of integers with each integer which is the order of PT in each SoC
  Tdat2['b']<-as.integer(factor(Tdat2$AEBODSYS))
  Tdat2['i']<-as.integer(factor(Tdat2$AEDECOD))
  Tdat2<-Tdat2 %>% group_by(b) %>% mutate(j = as.integer(factor(AEDECOD)))

  # make sure Tdat2 is order by b and j for further analysis
  Tdat2<-Tdat2[order(Tdat2$b, Tdat2$j), ]

  # return the dataframe
  return (Tdat2)
}
