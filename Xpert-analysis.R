#######################################
###### File: Xpert-analysis.R #########
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/25 #############
###### Updated 2024/05/07 #############
#######################################

require(readxl)
require(tidyverse)

source("A_Const.R")
source("Sigmas.R")
source("Full_Analysis.R")

## Read in (simulated) data from data folder:
load("data/Xpert-data-sim.Rda")

## Get unique periods, clusters, start times, N, J:
Periods <- unique(xpert.dat$Period)
J <- length(Periods)

StartTimes <- xpert.dat %>% dplyr::filter(Interv==1) %>%
  group_by(Cluster) %>% dplyr::summarize(StartPd=min(Period)) %>%
  dplyr::arrange(StartPd,Cluster)
N <- length(StartTimes$Cluster)

## Prep Outcome Data in appropriate order:
Ord_Data <- xpert.dat %>% left_join(StartTimes, by="Cluster") %>%
  dplyr::arrange(StartPd,Cluster,Period)
Obs_Y <- matrix(data=c(Ord_Data$Outcome, Ord_Data$logOdds), ncol=2)
colnames(Obs_Y) <- c("Probability","Log Odds")

## Generate A matrix:
Amat <- gen_A(N,J)

## Run Solver for different assumption settings:
SO2 <- Solve_Assumption(Amat,StartTimes,J,
                        Assumption=2,
                        v.Mat=cbind(Avg=c(rep(1/28,28)),
                                    AvgEx7=c(rep(1/21,21),rep(0,7)),
                                    D.1=c(1/6,1/6,0,1/6,0,0,1/6,0,0,0,1/6,0,0,0,0,1/6,rep(0,12)),
                                    D.2=c(0,0,1/5,0,1/5,0,0,1/5,0,0,0,1/5,0,0,0,0,1/5,rep(0,11)),
                                    D.12=c(1/11,1/11,1/11,1/11,1/11,0,1/11,1/11,0,0,1/11,1/11,0,0,0,1/11,1/11,rep(0,11)),
                                    T.234=c(rep(1/6,6),rep(0,22)),
                                    U.1=c(1,rep(0,27)),
                                    U.2=c(0,1,rep(0,26)),
                                    U.3=c(0,0,1,rep(0,25)),
                                    U.4=c(rep(0,3),1,rep(0,24)),
                                    U.5=c(rep(0,4),1,rep(0,23)),
                                    U.6=c(rep(0,5),1,rep(0,22)),
                                    U.7=c(rep(0,6),1,rep(0,21)),
                                    U.8=c(rep(0,7),1,rep(0,20)),
                                    U.9=c(rep(0,8),1,rep(0,19)),
                                    U.10=c(rep(0,9),1,rep(0,18)),
                                    U.11=c(rep(0,10),1,rep(0,17)),
                                    U.12=c(rep(0,11),1,rep(0,16)),
                                    U.13=c(rep(0,12),1,rep(0,15)),
                                    U.14=c(rep(0,13),1,rep(0,14)),
                                    U.15=c(rep(0,14),1,rep(0,13)),
                                    U.16=c(rep(0,15),1,rep(0,12)),
                                    U.17=c(rep(0,16),1,rep(0,11)),
                                    U.18=c(rep(0,17),1,rep(0,10)),
                                    U.19=c(rep(0,18),1,rep(0,9)),
                                    U.20=c(rep(0,19),1,rep(0,8)),
                                    U.21=c(rep(0,20),1,rep(0,7))),
                        save_loc="../int_large/",
                        save_prefix="xpert-solve-a_")
SO3 <- Solve_Assumption(Amat,StartTimes,J,
                        Assumption=3,
                        v.Mat=cbind(Avg=rep(1/7,7),
                             AvgEx7=c(rep(1/6,6),0),
                             D.1=c(1,rep(0,6)),
                             D.2=c(0,1,rep(0,5)),
                             D.3=c(0,0,1,rep(0,4)),
                             D.4=c(0,0,0,1,rep(0,3)),
                             D.5=c(rep(0,4),1,0,0),
                             D.6=c(rep(0,5),1,0),
                             D.7=c(rep(0,6),1),
                             Middle=c(0,1/3,1/3,1/3,0,0,0)),
                        save_loc="../int_large/",
                        save_prefix="xpert-solve-a_")
SO4 <- Solve_Assumption(Amat,StartTimes,J,
                        Assumption=4,
                        v.Mat=cbind(AvgEx7=c(rep(1/6,6),0),
                             T.1=c(1,rep(0,6)),
                             T.2=c(0,1,rep(0,5)),
                             T.3=c(0,0,1,rep(0,4)),
                             T.4=c(0,0,0,1,rep(0,3)),
                             T.5=c(rep(0,4),1,0,0),
                             T.6=c(rep(0,5),1,0),
                             T.7=c(rep(0,6),1),
                             Middle=c(0,1/3,1/3,1/3,0,0,0)),
                        save_loc="../int_large/",
                        save_prefix="xpert-solve-a_")
SO5 <- Solve_Assumption(Amat,StartTimes,J,
                        Assumption=5,
                        v.Mat=1,
                        save_loc="../int_large/",
                        save_prefix="xpert-solve-a_")

## Run variance minimizer for different settings:

### Independence:
for (i in 2:5) {
  assign(x=paste0("MVOut_",i,"_Ind"),
         value=MV_Assumption(SolveOut=get(paste0("SO",i)),
                             Assumption=i,
                             Sigma=create_Sigma_Ind(N=N,J=J),
                             SigmaName="Ind",
                             Observations=Obs_Y,
                             Permutations=1000,
                             save_loc="int/",
                             save_prefix="xpert-mv-a_"))
}

### Exchangeable (rho = 0.003 from Thompson et al. 2018):
for (i in 2:5) {
  assign(x=paste0("MVOut_",i,"_CS_0_003"),
         value=MV_Assumption(SolveOut=get(paste0("SO",i)),
                             Assumption=i,
                             Sigma=create_Sigma_CS(rho=0.003,N=N,J=J),
                             SigmaName="CS_0_003",
                             Observations=Obs_Y,
                             Permutations=1000,
                             save_loc="int/",
                             save_prefix="xpert-mv-a_"))
}

### AR(1) (rho = 0.012 gives average ICC within a cluster ~0.003):
for (i in 2:5) {
  assign(x=paste0("MVOut_",i,"_AR1_0_012"),
         value=MV_Assumption(SolveOut=get(paste0("SO",i)),
                             Assumption=i,
                             Sigma=create_Sigma_AR1(rho=0.012,N=N,J=J),
                             SigmaName="AR1_0_012",
                             Observations=Obs_Y,
                             Permutations=1000,
                             save_loc="int/",
                             save_prefix="xpert-mv-a_"))
}

## Import Results:
Assns <- 2:5
SigmaNames <- c("Ind","CS_0_003","AR1_0_012")

for (j in SigmaNames) {
  for (i in Assns) {
    load(file=paste0("int/xpert-mv-a_",i,"_",j,".Rda"))
  }
}

## Summarize Results:
for (j in SigmaNames) {
  print(paste0("Variance ",j))
  for (i in Assns) {
    print(paste0("Assumption ",i))
    print((get(paste0("MVOut_",i,"_",j))[["MV"]])[["Variance"]])
  }
}

for (i in Assns) {
  print(paste0("Assumption ",i))
  ProbEsts <- NULL
  OREsts <- NULL
  for (j in SigmaNames) {
    ProbEsts <- cbind(ProbEsts,(get(paste0("MVOut_",i,"_",j))[["Estimates"]])[,"Probability"])
    OREsts <- cbind(OREsts,exp((get(paste0("MVOut_",i,"_",j))[["Estimates"]])[,"Log Odds"]))
  }
  colnames(ProbEsts) <- SigmaNames
  colnames(OREsts) <- SigmaNames
  print(ProbEsts)
  print(OREsts)
}


### Inference:
for (i in Assns) {
  print(paste0("P-Values for Assumption ",i))
  for (j in SigmaNames) {
    print(j)
    print((get(paste0("MVOut_",i,"_",j)))[["P_Values"]])
  }
}
