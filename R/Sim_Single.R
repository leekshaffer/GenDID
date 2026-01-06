#######################################
###### File: Sim_Single.R #############
###### Lee Kennedy-Shaffer ############
#######################################

source("R/Analysis.R")
source("R/Ext_Comps.R")

library(dplyr)
library(tidyr)
library(tibble)
library(did) ## For CS
library(fixest) ## For SA
library(DIDmultiplegt) ## For CH
library(lme4) ## For mixed effects models and CLWP

## Load Setup Settings and Simulated Data

load("int/sim-setup.Rda")
load("int/sim-mvo-list.Rda")
load("int/sim-CI-objects.Rda")

# Analyze_One function (similar to Sim_Analyze.R version)

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

Analyze_One <- function(Scen,
                        SimNo,
                        NumPerms=0,
                        MVO_list,
                        Comps.Nest=NULL,
                        CI.GDID=NULL,
                        CI.SV.Mesh=NULL,
                        CI.Perc=0.95,
                        Comps=c("TW","CS","SA","CH","MEM","CPI","CPI.T","CPI.D","CPI.DT","CLWP","CLWPA"),
                        Comps_PermPs=c("TW","CS","SA","CH","MEM","CPI","CPI.T","CPI.D","CPI.DT","CLWP","CLWPA")) {
  ## Load Data
  load(paste0("int/sim-data-",Scen,".Rda"))
  # if (is.null(Params)) {
  #   Params <- Param_Set %>% dplyr::filter(Scenario==Scen)
  #   if (SimNo > (Params %>% pull(NumSims))) {
  #     stop(simpleError("SimNo is higher than the number of simulations available."))
  #   }
  # }

  ### Simulation-Specific Data:

  Data <- get(paste0("sim_data_",Scen))[[SimNo]]
  N <- length(unique(Data$Cluster))
  J <- length(unique(Data$Period))
  Obs_Y <- as.matrix(Data %>% dplyr::select(Y.ij.bar) %>%
                       dplyr::rename(Probability=Y.ij.bar) %>%
                       dplyr::mutate(`Log Odds`=if_else(Probability==0,
                                                        log((0.5/101)/(1-0.5/101)),
                                                        log(Probability/(1-Probability)))))

  ### GenDID Results from An_From_Obj:

  An.Out <- An_From_Obj(Obs_Y=Obs_Y,
                        MVO_list=MVO_list,
                        NumPerms=NumPerms,
                        Comps.Nest=Comps.Nest,
                        CI.GDID=CI.GDID,
                        CI.SV.Mesh=CI.SV.Mesh,
                        CI.Perc=CI.Perc)
  Results <- An.Out$Results
  P.Orders <- An.Out$Permutation.Orders

  ### Comparisons from Comp_Outs:

  Data.Long <- Data %>%
    dplyr::rename(Summ=Y.ij.bar) %>%
    dplyr::select(-Y.ij.sd) %>%
    tidyr::pivot_longer(cols=starts_with("Y.ij."),
                        names_to="Indiv", names_prefix="Y.ij.",
                        values_to="Outcome")

  Comp_Outs <- Ext_Comps(Data.Long=Data.Long,
                         SummOutName="Summ",
                         Comps=Comps,
                         Comps_PermPs=Comps_PermPs,
                         P.Orders=P.Orders,
                         NumPerms=0,
                         Results=Results,
                         CI.Perc=CI.Perc)

  return(list(Results=Results,
              Comparisons=Comp_Outs))
}

# Speed Testing:

## From test on Work Comp: 200 seconds for 100 permutations w/ mesh of 501
## From test on Work Comp: 8 minutes for 250 permutations w/ mesh of 501
st <- proc.time()
CompFull <- c("TW","CS","SA","CH","MEM","CPI","CPI.T","CPI.D","CPI.DT","CLWP","CLWPA")
A1 <- Analyze_One(Scen=1,
                  SimNo=1,
                  NumPerms=100,
                  MVO_list=MVO_list_full,
                  Comps.Nest=Comp_wts,
                  CI.GDID=CI.Tx.Obj_1,
                  CI.SV.Mesh=Single_Vals,
                  CI.Perc=0.95,
                  Comps=CompFull,
                  Comps_PermPs=CompFull)
proc.time() - st
