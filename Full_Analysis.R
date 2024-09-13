#######################################
###### File: Full_Analysis.R ##########
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/26 #############
###### Updated 2024/08/12 #############
#######################################

source("D_F_Const.R")
source("Rank_Analysis.R")
source("Var_Min.R")
source("Solver.R")

## Inputs to Solve_Assumption function:
### Amat is the A matrix for the setting
### StartTimes is a data frame/tibble with two columns that will be input to gen_DFT:
##### Cluster has the cluster names in order by starting time
##### StartPd has the starting times for those clusters
### OrderedPds is the vector of labels for the periods that will be observed, in temporal order
### Assumption is the numbered assumption setting for the treatment effects (1--5)
### v.Mat is the v vector of the desired estimand, or a matrix of such column vectors
### save_loc is the folder to save it in (be sure to end with / if anything other than blank is passed)
### save_prefix is the prefix to put on the save files

## Function to run solver for any assumption and v matrix
Solve_Assumption <- function(Amat,StartTimes,OrderedPds,
                             Assumption,v.Mat,
                             save_loc="",
                             save_prefix="solve-a_") {
  ## Generate D,F,Theta matrices:
  DFT_int <- gen_DFT(Clusters=StartTimes$Cluster,
                     StartPeriods=StartTimes$StartPd,
                     OrderedPds=OrderedPds,
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
       file=paste0(save_loc,save_prefix,Assumption,".Rda"))
  return(get(paste0("SolveOut_",Assumption)))
}

## Permutation Inference
## StartTimes is a data frame with a Cluster column and StartPd column
## J is the total number of observed periods
## data is a data frame with the outcomes, sorted by Cluster and then Period
Permute_data <- function(StartTimes, J, data, Obs.weights) {
  ST_P <- tibble(StartPd=rep(sample(StartTimes$StartPd, size=nrow(StartTimes), replace=FALSE), each=J),
                 Cluster=rep(sample(StartTimes$Cluster, size=nrow(StartTimes), replace=FALSE), each=J))
  return(t(Obs.weights) %*% as.matrix(cbind(ST_P, data) %>% dplyr::arrange(StartPd,Cluster) %>% 
                                        dplyr::select(-c(Cluster,StartPd))))
}
## Observations is a vector with the observed outcomes, ordered by cluster then period
## N is the total number of clusters/units
## J is the total number of periods (Note: N*J should be the number of rows in Observations)
## Obs.weights is the observation weights as output from min_var
Permute_order <- function(N, J) {
  return((rep(sample.int(N, size=N, replace=FALSE), each=J)-1)*J+rep(1:J, times=N))
}

Permute_obs <- function(Observations, N, J, Obs.weights) {
  Order <- Permute_order(N, J)
  if (is.null(Obs.weights)) {
    return(Observations[Order])
  } else {
    return(t(Obs.weights) %*% as.matrix(Observations)[Order,,drop=FALSE])
  }
}

## Inputs to MV_Assumption function:
### SolveOut is the output from a call to Solve_Assumption
### Assumption is the numbered assumption setting for the treatment effects (1--5)
### Sigma is the covariance matrix (up to scalar) to use for minimizing variance
### SigmaName (optional) is a short name to use in the files saved to indicate the variance setting
### Observations: if given, will compute estimates for the minimum-variance estimator using these observations
#### They must be in order of the A matrix:
#### by cluster in the order given in Solve_Assumption input StartTimes, then by period within cluster
### save_loc is the folder to save it in (be sure to end with / if anything other than blank is passed)
### save_prefix is the prefix to put on the save files
MV_Assumption <- function(SolveOut, Assumption,
                          Sigma, 
                          SigmaName=NULL,
                          Observations=NULL,
                          Permutations=NULL,
                          save_loc="",
                          save_prefix="mv-a_") {
  MV_int <- min_var(solve_obj=SolveOut$Solve,
                    A_mat=SolveOut$Amat,
                    Sigma=Sigma)
  if (is.null(Observations)) {
    D_Full <- cbind(SolveOut$DFT$D_aug,
                    MV_int$DID.weights)
    assign(x=paste0("MV_",Assumption,"_",SigmaName),
           value=list(Amat=SolveOut$Amat,
                            D_Full=D_Full,
                            MV=MV_int))
    save(list=paste0("MV_",Assumption,"_",SigmaName),
         file=paste0(save_loc,save_prefix,Assumption,"_",SigmaName,".Rda"))
    return(get(paste0("MV_",Assumption,"_",SigmaName)))
  } else {
    Estimates <- t(MV_int$Obs.weights) %*% Observations
    D_Full <- cbind(SolveOut$DFT$D_aug,
                    SolveOut$Amat %*% Observations,
                    MV_int$DID.weights)
    if(!is.null(Permutations)) {
      Perms <- replicate(n=Permutations,
                         expr=Permute_obs(Observations=Observations, 
                                          N=SolveOut$DFT$N, J=SolveOut$DFT$J, 
                                          Obs.weights=MV_int$Obs.weights))
      PermRes <- simplify2array(apply(Perms, 3, 
                                      FUN=function(x) abs(x) >= abs(Estimates), simplify=FALSE))
      PVals <- apply(PermRes, c(1,2), mean)
    } else {
      PVals <- NULL
    }
    assign(x=paste0("MVOut_",Assumption,"_",SigmaName),
           value=list(Amat=SolveOut$Amat,
                      Estimates=Estimates,
                      D_Full=D_Full,
                      MV=MV_int,
                      P_Values=PVals))
    save(list=paste0("MVOut_",Assumption,"_",SigmaName), 
         file=paste0(save_loc,save_prefix,Assumption,"_",SigmaName,".Rda"))
    return(get(paste0("MVOut_",Assumption,"_",SigmaName)))
  }
}



