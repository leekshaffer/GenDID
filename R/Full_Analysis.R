#######################################
###### File: Full_Analysis.R ##########
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/26 #############
###### Updated 2024/08/12 #############
#######################################

source("R/D_F_Const.R")
source("R/Rank_Analysis.R")
source("R/Var_Min.R")
source("R/Solver.R")

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
    return(Observations[Order,])
  } else {
    return(list(Obs=Observations[Order,],
                Ests=t(Obs.weights) %*% as.matrix(Observations)[Order,,drop=FALSE]))
  }
}

#' Compute Minimum-Variance Estimator and Save Results
#'
#' This function computes the minimum-variance estimator for treatment effects
#' under specified assumptions and optionally computes estimates using provided
#' observations. Results are saved to files for further analysis.
#'
#' @param SolveOut A list output from a call to `Solve_Assumption`, containing:
#'   - `Solve`: A solve object.
#'   - `Amat`: The A matrix used in the computation.
#'   - `DFT`: A list with elements `D_aug`, `N`, and `J` used for computations.
#' @param Assumption An integer (1–5) indicating the assumption setting for treatment effects.
#' @param Sigma The covariance matrix (up to scalar) to use for minimizing variance.
#' @param SigmaName A short string to indicate the variance setting in saved files (optional).
#' @param Observations A matrix of observations for computing estimates of the minimum-variance
#'   estimator (optional). Observations must be ordered by cluster and period as in the `StartTimes`
#'   input to `Solve_Assumption`.
#' @param Permutations An integer specifying the number of permutations to compute p-values
#'   for the observed estimates (optional).
#' @param save_loc A string specifying the folder path to save the output. Ensure the string ends
#'   with a "/" if it is not blank.
#' @param save_prefix A string prefix for the filenames of the saved results (default: "mv-a_").
#'
#' @return A list containing:
#'   - `Amat`: The A matrix used in the computation.
#'   - `D_Full`: A matrix combining augmented D values, optionally with observations and weights.
#'   - `MV`: The result of the `min_var` computation.
#'   - `Estimates` (if `Observations` is provided): A matrix of estimated treatment effects.
#'   - `P_Values` (if `Permutations` is provided): A matrix of permutation-based p-values.
#'
#' @details
#' If `Observations` is not provided, the function computes the minimum-variance estimator
#' weights and saves the results. If `Observations` are provided, the function computes
#' estimates for the minimum-variance estimator and optionally performs permutation testing
#' to calculate p-values.
#'
#' Saved files are named using the format:
#' `<save_prefix><Assumption>_<SigmaName>.Rda`
#'
#' @examples
#' # Example usage:
#' result <- MV_Assumption(SolveOut = solve_out,
#'                         Assumption = 3,
#'                         Sigma = diag(10),
#'                         SigmaName = "example",
#'                         Observations = obs_matrix,
#'                         Permutations = 100,
#'                         save_loc = "results/",
#'                         save_prefix = "mv-a_")
#'
#' @export
MV_Assumption <- function(SolveOut,
                          Assumption,
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



