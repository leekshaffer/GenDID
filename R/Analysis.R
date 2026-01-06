#######################################
###### File: Analysis.R ###############
###### Lee Kennedy-Shaffer ############
#######################################

source("R/Var_Min.R")

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

# An_From_Obj helper function

### Inputs:
#### Scen: simulation scenario number
#### SimNo: the number of the simulation within that scenario
#### NumPerms: the number of permutations to do for P-values and CIs
#### MVO_list: List of min_var objects to use for analysis; must be a
####  named list where the names will be used as prefixes for the estimators
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
    Obs.W <- cbind(Obs.W, Comp_wts)
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





