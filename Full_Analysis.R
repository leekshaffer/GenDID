#######################################
###### File: Full_Analysis.R ##########
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/26 #############
###### Updated 2024/04/26 #############
#######################################

source("D_F_Const.R")
source("Rank_Analysis.R")
source("Var_Min.R")
source("Solver.R")

## Read in data:
load("data/Xpert-data.Rda")

## Inputs to Solve_Assumption function:
### Amat is the A matrix for the setting
### StartTimes is a data frame/tibble with two columns that will be input to gen_DFT:
##### Cluster has the cluster names in order by starting time
##### StartPd has the starting times for those clusters
### J is the number of total periods
### Assumption is the numbered assumption setting for the treatment effects (1--5)
### v.Mat is the v vector of the desired estimand, or a matrix of such column vectors

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

## Inputs to MV_Assumption function:
### SolveOut is the output from a call to Solve_Assumption
### Assumption is the numbered assumption setting for the treatment effects (1--5)
### Sigma is the covariance matrix (up to scalar) to use for minimizing variance
### SigmaName (optional) is a short name to use in the files saved to indicate the variance setting
### Observations: if given, will compute estimates for the minimum-variance estimator using these observations
#### They must be in order of the A matrix:
#### by cluster in the order given in Solve_Assumption input StartTimes, then by period within cluster
MV_Assumption <- function(SolveOut, Assumption,
                          Sigma, 
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