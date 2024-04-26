#######################################
###### File: Xpert-analysis.R #########
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/25 #############
###### Updated 2024/04/25 #############
#######################################

require(readxl)
require(tidyverse)

source("A_Const.R")
source("D_F_Const.R")
source("Sigmas.R")
source("Rank_Analysis.R")
source("Var_Min.R")
source("Solver.R")

## Read in data:
load("data/Xpert-data.Rda")

## Get unique periods, clusters, start times, N, J:
Periods <- unique(trajclust2$Period)
J <- length(Periods)

StartTimes <- trajclust2 %>% dplyr::filter(Interv==1) %>%
  group_by(Cluster) %>% dplyr::summarize(StartPd=min(Period)) %>%
  dplyr::arrange(StartPd,Cluster)
N <- length(StartTimes$Cluster)

## Prep Outcome Data in appropriate order:
Ord_Data <- trajclust2 %>% left_join(StartTimes, by="Cluster") %>%
  dplyr::arrange(StartPd,Cluster,Period)
Obs_Y <- matrix(data=c(Ord_Data$Outcome, Ord_Data$logOdds), ncol=2)
colnames(Obs_Y) <- c("Probability","Log Odds")

## Generate A matrix:
Amat <- gen_A(N,J)

## Function to run solver for any assumption and v matrix
Solve_Assumption <- function(Amat,StartTimes,J,
                             Assumption,v.Mat) {
  ## Generate D,F,Theta matrices:
  DFT_int <- gen_DFT(Clusters=StartTimes$Cluster,
                     StartPeriods=StartTimes$StartPd,
                     J=J,
                     Assumption=Assumption)
  rank_int <- rank_an(DFT_obj=DFT_int,
                      v=v.Mat)
  solve_int <- solve_WA(DFT_obj=DFT_int,
                        A_mat=Amat,
                        v=v.Mat,
                        rank_obj=rank_int,
                        DID_full=TRUE)
  assign(x=paste0("SolveOut_",Assumption),
         value=list(Amat=Amat,
                    DFT=DFT_int,
                    Solve=solve_int))
  save(list=paste0("SolveOut_",Assumption), 
       file=paste0("../int_large/xpert-solve-a_",Assumption,".Rda"))
  return(get(paste0("SolveOut_",Assumption)))
}

## Function to run var minimizer for any assumption and v matrix
MV_Assumption <- function(SolveOut, Sigma, Assumption,
                          SigmaName=NULL,
                          Observations=NULL) {
  MV_int <- min_var(solve_obj=SolveOut$Solve,
                    A_mat=SolveOut$Amat,
                    Sigma=Sigma)
  if (is.null(Observations)) {
    assign(x=paste0("MVOut_",Assumption,"_",SigmaName),
           value=MV_int)
    save(list=paste0("MV_",Assumption,"_",SigmaName),
         file=paste0("../int_large/xpert-mv-a_",Assumption,"_",SigmaName,".Rda"))
    return(MV_int)
  } else {
    Estimates <- t(MV_int$Obs.weights) %*% Observations
    D_Full <- cbind(SolveOut$DFT$D_aug,
                    SolveOut$Amat %*% Observations,
                    MV_int$DID.weights)
    assign(x=paste0("MVOut_",Assumption,"_",SigmaName),
           value=list(Amat=Amat,
                      Estimates=Estimates,
                      D_Full=D_Full,
                      MV=MV_int))
    save(list=paste0("MVOut_",Assumption,"_",SigmaName), 
         file=paste0("int/xpert-mv-a_",Assumption,"_",SigmaName,".Rda"))
    return(get(paste0("MVOut_",Assumption,"_",SigmaName)))
  }
}

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
                             Middle=c(0,1/3,1/3,1/3,0,0,0)))
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
                             Middle=c(0,1/3,1/3,1/3,0,0,0)))
SO5 <- Solve_Assumption(Amat,StartTimes,J,
                        Assumption=5,
                        v.Mat=1)
## Timing (desktop) notes: <1min each for A3/A4/A5

MV3 <- MV_Assumption(SolveOut=SO3,
                     Sigma=create_Sigma_CS(rho=0.003,N=N,J=J),
                     Assumption=3,
                     SigmaName="CS0_003",
                     Observations=Obs_Y)
MV4 <- MV_Assumption(SolveOut=SO4,
                     Sigma=create_Sigma_CS(rho=0.003,N=N,J=J),
                     Assumption=4,
                     SigmaName="CS0_003",
                     Observations=Obs_Y)
MV5 <- MV_Assumption(SolveOut=SO5,
                     Sigma=create_Sigma_CS(rho=0.003,N=N,J=J),
                     Assumption=5,
                     SigmaName="CS0_003",
                     Observations=Obs_Y)
## Timing (desktop) notes: ~3sec each for A3/A4/A5
