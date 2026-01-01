#######################################
###### File: Full_Analysis.R ##########
###### Lee Kennedy-Shaffer ############
#######################################

require(dplyr)
require(tibble)
require(tidyr)
require(stringr)

# source("R/Mat_Const.R")
# source("R/Rank_Analysis.R")
# source("R/Var_Min.R")
# source("R/Solver.R")
# source("R/Permuation_Fns.R")

# MV_Assumption function

### Inputs:
#### SolveOut: The output from Solve_Assumption.
#### Assumption: Number 1--5 corresponding to the Assumption Settings.
#### Sigma: The variance-covariance matrix to use for minimizing variance.
#### SigmaName: A short string to indicate the variance setting in saved files (optional).
#### Observations: A matrix of observations for computing estimates of the estimator;
####  Ordered by unit and then period.
#### Permutations: An integer specifying the number of permutations to compute p-values
####  for the observed estimates (optional).
### Output: A list of the following:
#### A_mat: The A matrix used in the computation.
#### Estimates: A matrix of estimated treatment effects (if Observations is provided)
#### D_Full: A matrix combining augmented D values, optionally with observations and weights.
#### MV: The output of the min_var function.
#### P_Values: A matrix of permutation-based p-values (if Permutations is provided)

MV_Assumption <- function(SolveOut,
                          Assumption,
                          Sigma,
                          SigmaName=NULL,
                          Observations=NULL,
                          Permutations=NULL,
                          save_loc="",
                          save_prefix="mv-a_") {
  MV_int <- min_var(solve_obj=SolveOut$Solve,
                    A_mat=SolveOut$ADFT$A_mat,
                    Sigma=Sigma)
  if (is.null(Observations)) {
    D_Full <- cbind(SolveOut$DFT$D_aug,
                    MV_int$DID.weights)
    assign(x=paste0("MV_",Assumption,"_",SigmaName),
           value=list(Amat=SolveOut$ADFT$A_mat,
                      D_Full=D_Full,
                      MV=MV_int))
    # TODO: Inform user that you are saving variable
    # to a .Rda file using the specified save_loc and save_prefix.
    save(list=paste0("MV_",Assumption,"_",SigmaName),
         file=paste0(save_loc,save_prefix,Assumption,"_",SigmaName,".Rda"))
    return(get(paste0("MV_",Assumption,"_",SigmaName)))
  } else {
    Estimates <- t(MV_int$Obs.weights) %*% Observations
    D_Full <- cbind(SolveOut$ADFT$D_aug,
                    SolveOut$ADFT$A_mat %*% Observations,
                    MV_int$DID.weights)
    if(!is.null(Permutations)) {
      Perms <- replicate(n=Permutations,
                         expr=Permute_obs(Observations=Observations,
                                          N=SolveOut$ADFT$N, J=SolveOut$ADFT$J,
                                          Obs.weights=MV_int$Obs.weights))
      PermRes <- simplify2array(apply(Perms, 3,
                                      FUN=function(x) abs(x) >= abs(Estimates), simplify=FALSE))
      PVals <- apply(PermRes, c(1,2), mean)
    } else {
      PVals <- NULL
    }
    assign(x=paste0("MVOut_",Assumption,"_",SigmaName),
           value=list(A_mat=SolveOut$ADFT$A_mat,
                      Estimates=Estimates,
                      D_Full=D_Full,
                      MV=MV_int,
                      P_Values=PVals))
    # TODO: Inform user that you are saving variable
    # to a .Rda file using the specified save_loc and save_prefix.
    # save(list=paste0("MVOut_",Assumption,"_",SigmaName),
    #      file=paste0(save_loc,save_prefix,Assumption,"_",SigmaName,".Rda"))
    return(get(paste0("MVOut_",Assumption,"_",SigmaName)))
  }
}
