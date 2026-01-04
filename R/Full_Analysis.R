#######################################
###### File: Full_Analysis.R ##########
###### Lee Kennedy-Shaffer ############
#######################################

library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

source("R/Mat_Const.R")
source("R/Sigmas.R")
source("R/Rank_Analysis.R")
source("R/Var_Min.R")
source("R/Solver.R")
source("R/Permutation_Fns.R")

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
#### CI_list: a (preferably named) list of vector or matrix of column vectors of
####  treatment effects to test for inclusion in CI (optional). Requires Permutations as well.
### Output: A list of the following:
#### A_mat: The A matrix used in the computation (if Observations is provided)
#### ADFT: The ADFT object from gen_ADFT (if Observations is not provided)
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
                          CI_list=NULL #,
                          # save_loc="",
                          # save_prefix="mv-a_"
                          ) {
  MV_int <- min_var(solve_obj=SolveOut$Solve,
                    A_mat=SolveOut$ADFT$A_mat,
                    Sigma=Sigma)
  if (is.null(Observations)) {
    D_Full <- cbind(SolveOut$ADFT$D_aug,
                    MV_int$DID.weights)
    assign(x=paste0("MV_",Assumption,"_",SigmaName),
           value=list(ADFT=SolveOut$ADFT,
                      D_Full=D_Full,
                      MV=MV_int))
    return(get(paste0("MV_",Assumption,"_",SigmaName)))
  } else {
    Estimates <- t(MV_int$Obs.weights) %*% Observations
    D_Full <- cbind(SolveOut$ADFT$D_aug,
                    SolveOut$ADFT$A_mat %*% Observations,
                    MV_int$DID.weights)
    if(!is.null(Permutations)) {
      if (is.null(CI_list)) {
        Perms <- replicate(n=Permutations,
                           expr=Permute_obs(Observations=Observations,
                                            N=SolveOut$ADFT$N, J=SolveOut$ADFT$J,
                                            Obs.weights=MV_int$Obs.weights)$Ests)
        PermRes <- simplify2array(apply(Perms, 3,
                                        FUN=function(x) abs(x) >= abs(Estimates), simplify=FALSE))
        PVals <- apply(PermRes, c(1,2), mean)
        CI_Checks <- NULL
      } else {
        CI_Tx_All <- CI_get_Tx_All(SolveOut$ADFT, CI_list)
        CI_Ests <- c(list(Estimates=Estimates),
                     lapply(CI_Tx_All, function(x) Estimates - t(MV_int$Obs.weights) %*% as.matrix(x)))
        Perms <- replicate(n=Permutations,
                           expr=Permute_obs(Observations=Observations,
                                            N=SolveOut$ADFT$N, J=SolveOut$ADFT$J,
                                            Obs.weights=MV_int$Obs.weights,
                                            CI.Tx.List=CI_Tx_All,
                                            Drop_Obs=TRUE))
        Outlist <- NULL
        for (i in 1:length(CI_Ests)) {
          PermsRow <- Perms[i,]
          PermsRowRes <- simplify2array(lapply(PermsRow, function(x) abs(x) >= abs(CI_Ests[[i]])))
          Outlist[[i]] <- apply(PermsRowRes, c(1,2), mean)
        }
        names(Outlist) <- names(CI_Ests)
        PVals <- Outlist[["Estimates"]]
        CI_Checks <- Outlist[names(CI_Tx_All)]
      }
    } else {
      PVals <- NULL
      CI_Checks <- NULL
    }
    assign(x=paste0("MVOut_",Assumption,"_",SigmaName),
           value=c(list(A_mat=SolveOut$ADFT$A_mat,
                      Estimates=Estimates,
                      D_Full=D_Full,
                      MV=MV_int,
                      P_Values=PVals),
                   CI_Checks))
    # TODO: Inform user that you are saving variable
    # to a .Rda file using the specified save_loc and save_prefix.
    # save(list=paste0("MVOut_",Assumption,"_",SigmaName),
    #      file=paste0(save_loc,save_prefix,Assumption,"_",SigmaName,".Rda"))
    return(get(paste0("MVOut_",Assumption,"_",SigmaName)))
  }
}
