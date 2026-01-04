#######################################
###### File: Sim_Analyze.R ############
###### Lee Kennedy-Shaffer ############
#######################################

source("R/Full_Analysis.R")

set.seed(912350)

## Load Setup Settings and Simulated Data

load("int/sim-setup.Rda")
load("int/sim-mvo-list.Rda")
load("int/sim-CI-objects.Rda")
for (i in 1:(dim(Param_Set)[1])) {
  load(paste0("int/sim-data-",i,".Rda"))
}

# Analyze_One function

### Inputs:
#### Scenario: simulation scenario number
#### Params: the parameter set for the scenario; can be left as NULL
#### SimNo: the number of the simulation within that scenario
#### MVO_list: List of min_var objects to use for analysis
### Outputs: List of the following:
#### Estimates: Estimated treatment effects from the scenario/simulation
#### P_Values: P-Values for the effects from the scenario/simulation
#### CI_Ps: P-Values for perm. CIs from the scenario/simulation

Analyze_One <- function(Scen,
                        Params=NULL,
                        SimNo,
                        MVO_list) {
  if (is.null(Params)) {
    Params <- Param_Set %>% dplyr::filter(Scenario==Scen)
    if (SimNo > (Params %>% pull(NumSims))) {
      stop(simpleError("SimNo is higher than the number of simulations available."))
    }
  }
  CI.T.O <- get(paste0("CI.Tx.Obj_",Scen))
  Data <- get(paste0("sim_data_",Scen))[[SimNo]]
  Obs_Y <- as.matrix(Data %>% dplyr::select(Y.ij.bar) %>%
    dplyr::rename(Probability=Y.ij.bar) %>%
    dplyr::mutate(`Log Odds`=if_else(Probability==0,
                                     log((0.5/100)/(1-0.5/100)),
                                     log(Probability/(1-Probability)))))

  Estimates <- NULL
  Outlist <- NULL
  for (Setting in names(MVO_list)) {
    Obs.W <- MVO_list[[Setting]]$MV$Obs.weights
    if (is.null(colnames(Obs.W))) {
      colnames(Obs.W) <- Setting
    } else {
      colnames(Obs.W) <- paste(Setting, colnames(Obs.W), sep="_")
    }

    Ests <- t(Obs.W) %*% Obs_Y
    Estimates <- rbind(Estimates, Ests)

    CI_Ests <- c(list(PValues=Ests),
                 lapply(CI.T.O, function(x) Ests - t(Obs.W) %*% as.matrix(x)))

    Perms <- replicate(n=Params %>% pull(NumPerms),
                       expr=Permute_obs(Observations=Obs_Y,
                                        N=Params %>% pull(N),
                                        J=Params %>% pull(J),
                                        Obs.weights=Obs.W,
                                        CI.Tx.List=CI.T.O,
                                        Drop_Obs=TRUE))
    for (i in 1:length(CI_Ests)) {
      PermsRow <- Perms[i,]
      PermsRowRes <- simplify2array(lapply(PermsRow, function(x) abs(x) >= abs(CI_Ests[[i]])))
      if (length(Outlist) < i) {
        Outlist <- c(Outlist, list(apply(PermsRowRes, c(1,2), mean)))
      } else {
        Outlist[[i]] <- rbind(Outlist[[i]],
                              apply(PermsRowRes, c(1,2), mean))
      }
    }
  }

  names(Outlist) <- names(CI_Ests)
  # PVals <- Outlist[["Estimates"]]
  # CI_Checks <- Outlist[names(CI.T.O)]

  return(c(list(Estimates=Estimates),
           Outlist))
}

# Analyze_Scen function

### Inputs:
#### Scenario: simulation scenario number
#### MVO_list: List of min_var objects to use for analysis
### Outputs: Data frame of the summaries of estimates, p-values, and CIs from the simulations for that scenario

Analyze_Scen <- function(Scen,
                         MVO_list) {
  Params <- Param_Set %>% dplyr::filter(Scenario==Scen)
  List <- sapply(1:(Params %>% pull(NumSims)),
                 FUN=function(x) simplify2array(Analyze_One(Scen, Params,
                                             x, MVO_list)),
                 simplify="array")
  CI_names <- dimnames(List)[[3]][3:(dim(List)[3])]

  MeanEst <- apply(X=List[,,"Estimates",],
                   MARGIN=c(1,2),
                   FUN=mean)
  MedianEst <- apply(X=List[,,"Estimates",],
                     MARGIN=c(1,2),
                     FUN=median)
  SDEst <- apply(X=List[,,"Estimates",],
                     MARGIN=c(1,2),
                     FUN=sd)
  Power <- apply(X=List[,,"PValues",],
                 MARGIN=c(1,2),
                 FUN=function(x) mean(x < 0.05))
  Values <- rbind(t(MeanEst), t(MedianEst), t(SDEst), t(Power))
  for (name in CI_names) {
    Values <- rbind(Values,
                    t(apply(X=List[,,name,],
                          MARGIN=c(1,2),
                          FUN=function(x) mean(x >= 0.05))))
  }
  rownames(Values) <- NULL

  return(tibble(Scenario=rep(Scen, dim(List)[2]*(dim(List)[3]+2)),
                           Result=rep(c("Mean Estimate", "Median Estimate", "SD Estimate",
                                        "Power", paste("CI Coverage", CI_names)),
                                      each=dim(List)[2]),
                           Outcome=rep(dimnames(List)[[2]], dim(List)[3]+2)) %>%
    cbind(Values))
}

## Running Analysis:

### Non-Parallelized Version:

# Full_Sim_CI_Res <- do.call("rbind",
#                            lapply(1:9,
#                                   FUN=function(x) Analyze_Scen(Scen=x,
#                                                                MVO_list=MVO_list_full)))

### Parallelized Version:

library(foreach)
library(doParallel)

par_clust <- makeCluster(5)
registerDoParallel(par_clust)
Full_Sim_CI_Res <-
  foreach (x=Param_Set$Scenario,
           .combine=rbind,
           .inorder=TRUE) %dopar% {
    source("R/Full_Analysis.R")
    load("int/sim-setup.Rda")
    load("int/sim-mvo-list.Rda")
    load("int/sim-CI-objects.Rda")
    load(paste0("int/sim-data-",x,".Rda"))
    Analyze_Scen(Scen=x,
                 MVO_list=MVO_list_full)
  }
stopCluster(cl = par_clust)

save(Full_Sim_CI_Res,
     file=paste0("res","/Full_Sim_CI_Res"))
