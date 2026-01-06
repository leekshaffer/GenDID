#######################################
###### File: Analysis.R ###############
###### Lee Kennedy-Shaffer ############
#######################################

source("R/Mat_Const.R")
source("R/Solver.R")
source("R/Var_Min.R")
source("R/CI_Prep_Fns.R")

require(tibble)
require(tidyr)
require(dplyr)

# Permute_order helper function

### Inputs:
#### N: total number of units
#### J: total number of observed periods
### Output: randomly permuted unit order

Permute_order <- function(N, J) {
  return((rep(sample.int(N, size=N, replace=FALSE), each=J)-1)*J+rep(1:J, times=N))
}

# An_From_Obj function

### Inputs:
#### Obs_Y: the matrix of observation vectors to use
#### MVO_list: List of min_var objects to use for analysis; must be a
####  named list where the names will be used as prefixes for the estimators
#### NumPerms: the number of permutations to do for P-values and CIs
#### Comps.Nest: the data frame of observation weights for computation methods
####  that can be nested in GenDID framework
####  (e.g., $Obs.weights of Comp_Ests_Weights function); or NULL
#### CI.GDID: The CI effects object with effects to test
####  (e.g., the CI_get_Tx_All function output); or NULL
#### CI.SV.Mesh: The named list of the mesh of CI values,
####  corresponding to elements of CI.GDID; or NULL
#### CI.Perc: quantile to use for CI (e.g., 0.95 for 95% CI)
### Outputs: List of the following:
#### Estimates: Estimated treatment effects from the scenario/simulation
#### P_Values: P-Values for the effects from the scenario/simulation
#### CI_Ps: P-Values for perm. CIs from the scenario/simulation

An_From_Obj <- function(Obs_Y,
                        MVO_list,
                        NumPerms=0,
                        Comps.Nest=NULL,
                        CI.GDID=NULL,
                        CI.SV.Mesh=NULL,
                        CI.Perc=0.95) {
  N <- MVO_list[[1]]$ADFT$N
  J <- MVO_list[[1]]$ADFT$J

  ### Permuted Orders:
  if (NumPerms > 0) {
    P.Orders <- replicate(n=NumPerms,
                          expr=Permute_order(N, J),
                          simplify=FALSE)
  } else {
    P.Orders <- NULL
  }

  ## Get GenDID and Nested Method Estimates and Inference

  ### Full list of observation weights for assumption settings
  Obs.W <- do.call("cbind",
                   lapply(1:length(MVO_list),
                          function(x) as_tibble(MVO_list[[x]]$MV$Obs.weights) %>%
                            dplyr::rename_with(~paste(names(MVO_list)[x],
                                                      colnames(MVO_list[[x]]$MV$Obs.weights),
                                                      sep="_",recycle0=FALSE))))
  if (!is.null(Comps.Nest)) {
    Obs.W <- cbind(Obs.W, Comps.Nest)
  }

  ### Estimates calculation:
  Estimates <- t(Obs.W) %*% Obs_Y
  Results <- as_tibble(Estimates, rownames="Estimator") %>%
    tidyr::pivot_longer(cols=all_of(colnames(Obs_Y)),
                        names_to="Outcome",
                        values_to="Estimate")

  ### P Values calculation:
  if (!is.null(P.Orders)) {
    Perm_Ests <- lapply(P.Orders, FUN=function(x) t(Obs.W) %*% Obs_Y[x,])
    PVals <- apply(simplify2array(
      lapply(Perm_Ests,
             FUN=function(x) abs(x) >= abs(Estimates))),
      MARGIN=c(1,2), mean)
    Results <- Results %>%
      left_join(as_tibble(PVals, rownames="Estimator") %>%
                  tidyr::pivot_longer(cols=all_of(colnames(Obs_Y)),
                                      names_to="Outcome",
                                      values_to="P"),
                by=join_by(Estimator,Outcome))
  }

  ### CIs calculation:
  if (!is.null(P.Orders) & !is.null(CI.GDID)) {
    CI_Ps <- lapply(names(CI.GDID),
                    FUN=function(CI.name) {
                      CI.FX <- as.matrix(CI.GDID[[CI.name]])
                      CI.Estimates <- Estimates - t(Obs.W) %*% CI.FX
                      CI_Perm_Ests <- lapply(P.Orders,
                                             FUN=function(x) t(Obs.W) %*% (Obs_Y - CI.FX)[x,])
                      CI_Ps <- apply(simplify2array(
                        lapply(CI_Perm_Ests,
                               FUN=function(x) abs(x) >= abs(CI.Estimates))),
                        MARGIN=c(1,2), mean)
                    })
    names(CI_Ps) <- names(CI.GDID)
    if (is.null(CI.SV.Mesh)) {
      CI_List <- lapply(CI_Ps, FUN=function(x) x >= 1-CI.Perc)
    } else {
      CI_List <- lapply(CI_Ps[names(CI_Ps)[!(names(CI_Ps) %in% CI.SV.Mesh$Name)]],
                        FUN=function(x) x >= 1-CI.Perc)
      Results <- Results %>%
        left_join(do.call("rbind",
                          lapply(CI.SV.Mesh$Name,
                                 FUN=function(x) as_tibble(CI_Ps[[x]], rownames="Estimator") %>%
                                   dplyr::mutate(Name=x))) %>%
                    pivot_longer(cols=colnames(Obs_Y),
                                 names_to="Outcome", values_to="CI_P") %>%
                    left_join(CI.SV.Mesh %>% pivot_longer(cols=colnames(Obs_Y),
                                                          names_to="Outcome", values_to="Effect"),
                              by=join_by(Name,Outcome)) %>%
                    dplyr::group_by(Estimator,Outcome) %>%
                    dplyr::filter(CI_P >= 1-CI.Perc) %>%
                    dplyr::summarize(CIL=min(Effect),
                                     CIU=max(Effect),
                                     .groups="drop"),
                  by=join_by(Estimator,Outcome)) %>%
        dplyr::mutate(CIW=CIU-CIL)
    }
    for (x in names(CI_List)) {
      Results <- Results %>%
        left_join(as_tibble(CI_List[[x]], rownames="Estimator") %>%
                    tidyr::pivot_longer(cols=colnames(Obs_Y),
                                        names_to="Outcome",
                                        values_to=paste0("CI_",x)),
                  by=join_by(Estimator, Outcome))
    }
  }
  return(list(Results=Results,
              Permutation.Orders=P.Orders))
}

# An_Full function

### Inputs:
#### Data: No specific order required; needs Cluster, Period, Interv columns and the columns named in Outcome_Cols
#### OrderedPds: A vector of the Period names/numbers, in their temporal order from earliest to latest
#### Outcome_Cols: The names of columns to use as outcomes
#### Assumption: The assumption setting, an integer from 2 to 5
#### v.Mat: A matrix of treatment effect vectors (as columns) to use. Include colnames to get named output
#### SigmaList: A named list of variance matrix assumptions to use (e.g., output from create_Sigma_* functions).
#### Nested_Comps: a character vector of weight-based comparators
####  (options: "TW","CS","SA","CH","CO","NP")
#### NumPerms: The number of permutations to use for inference; set to 0 for no inference
#### CIs_List: A named list of specific CI values to check. Each element of the list
####  should be a tibble with columns corresponding to Outcome_Cols and a number of rows corresponding
####  to the unique treatmente effects in Theta (if fewer/more, will repeat/truncate with a warning)
#### CI.SV.Mesh: A tibble with columns Name and the Outcome_Cols with the values to check as a mesh.
####  The tibble must be ordered from smallest to largest value checked (with all columns in order).
#### CI.Perc: quantile to use for CI (e.g., 0.95 for 95% CI)
#### Keep_MVO: logial indicating whether the MVO list should be exported
### Outputs: List of the following:
#### Results: Tibble of estimates, p-values, and CIs for all variance settings and estimators
#### P.Orders: List of order vectors used in permutation (for reproducibility with external)
#### MVO: MVO list object (if Keep_MVO==TRUE)

An_Full <- function(Data,
                    OrderedPds,
                    Outcome_Cols,
                    Assumption,
                    v.Mat,
                    SigmaList,
                    Nested_Comps,
                    NumPerms=0,
                    CIs_List=NULL,
                    CI.SV.Mesh=NULL,
                    CI.Perc=0.95,
                    Keep_MVO=FALSE
) {
  Pd.Nums <- tibble(Period=OrderedPds, Num=1:length(OrderedPds))
  STs <- Data %>% dplyr::left_join(Pd.Nums, by=join_by(Period)) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize(StartPd=min(Num[Interv==1]))
  Ord_Data <- Data %>% left_join(STs, by=join_by(Cluster)) %>%
    arrange(StartPd,Cluster,Period)
  ADFT <- gen_ADFT(Clusters=STs$Cluster,
                   StartPeriods=STs$StartPd,
                   OrderedPds=OrderedPds,
                   Assumption=Assumption)
  SO <- Solve_Assumption(StartTimes=STs %>% dplyr::arrange(StartPd,Cluster),
                          OrderedPds=OrderedPds,
                          Assumption=Assumption,
                          v.Mat=v.Mat)
  MVO <- lapply(SigmaList,
                FUN=function(x) min_var(SolveOut=SO,
                                        Sigma=x,
                                        Drop_ADFT=FALSE,
                                        method="CG", maxit=10000))
  if (!is.null(Nested_Comps)) {
    Comp_Wts <- Comp_Ests_Weights(ADFT_obj=SO$ADFT,
                                  estimator=Nested_Comps)$Obs.weights
  } else {
    Comp_Wts <- NULL
  }

  if (!is.null(CI.SV.Mesh)) {
    SV_list <- split(CI.SV.Mesh %>% dplyr::select(all_of(Outcome_Cols)),
                     seq(nrow(CI.SV.Mesh)))
    names(SV_list) <- CI.SV.Mesh$Name
    CIs_List <- c(CIs_List, SV_list)
  }
  if (!is.null(CIs_List)) {
    CI.GDID <- CI_get_Tx_All(ADFT_obj=SO$ADFT,
                             CI_list=CIs_List)
  } else {
    CI.GDID <- NULL
  }

  if (Keep_MVO) {
    return(c(An_From_Obj(Obs_Y=as.matrix(Ord_Data %>% dplyr::select(all_of(Outcome_Cols))),
                       MVO_list=MVO,
                       NumPerms=NumPerms,
                       Comps.Nest=Comp_Wts,
                       CI.GDID=CI.GDID,
                       CI.SV.Mesh=CI.SV.Mesh,
                       CI.Perc=CI.Perc),
             list(MVO=MVO)))
  } else {
    return(An_From_Obj(Obs_Y=as.matrix(Ord_Data %>% dplyr::select(all_of(Outcome_Cols))),
                       MVO_list=MVO,
                       NumPerms=NumPerms,
                       Comps.Nest=Comp_Wts,
                       CI.GDID=CI.GDID,
                       CI.SV.Mesh=CI.SV.Mesh,
                       CI.Perc=CI.Perc))
  }
}




