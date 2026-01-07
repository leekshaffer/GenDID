#######################################
###### File: Sim_Cluster.R ############
#######################################

source("R/Analysis.R")
source("R/Ext_Comps.R")

ArrNo <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
JobName <- Sys.getenv('SLURM_JOB_NAME')

NumPerArr <- 10

library(dplyr)
library(tidyr)
library(tibble)
library(did) ## For CS
library(fixest) ## For SA
# library(DIDmultiplegt) ## For CH
library(lme4) ## For mixed effects models and CLWP
library(bbmle) ## For CLWP

## Load Setup Settings and Simulated Data

load("sim_data/sim-setup.Rda")
load("sim_data/sim-mvo-list.Rda")
load("sim_data/sim-CI-objects.Rda")

# Analyze_One function

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
  load(paste0("sim_data/sim_data_",Scen,".Rda"))

  ### Simulation-Specific Data:

  Data <- get(paste0("sim_data_",Scen))[[SimNo]]
  N <- length(unique(Data$Cluster))
  J <- length(unique(Data$Period))
  Obs_Y <- as.matrix(Data %>% dplyr::select(Summ) %>%
                       dplyr::rename(Probability=Summ) %>%
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

  ### Comparisons from Comp_Outs:

  Data.Long <- Data %>%
    dplyr::select(Cluster,Period,Start,Interv,Summ,starts_with("Y.ij.")) %>%
    tidyr::pivot_longer(cols=starts_with("Y.ij."),
                        names_to="Indiv", names_prefix="Y.ij.",
                        values_to="Outcome")

  Comp_Outs <- Ext_Comps(Data.Long=Data.Long,
                         SummOutName="Summ",
                         Comps=Comps,
                         Comps_PermPs=Comps_PermPs,
                         P.Orders=An.Out$Permutation.Orders,
                         NumPerms=NumPerms,
                         Results=An.Out$Results,
                         CI.Perc=CI.Perc)

  return(bind_rows(An.Out$Results,
                   Comp_Outs %>% dplyr::mutate(Outcome="Probability") %>%
                     dplyr::rename(Estimator=Method,
                                   P.Asy=P,
                                   P=P.Perm)))
}

m <- as.integer(substr(JobName, nchar(JobName), nchar(JobName)))
Sims <- ((ArrNo-1)*NumPerArr+1):(ArrNo*NumPerArr)
Perms_Use <- Param_Set$NumPerms[Param_Set$Scenario==m]

print(paste0("Scenario: ",m))
print(Sims)

Output <- NULL
for (SimVal in Sims) {
  print(paste0("SimVal: ",SimVal))
  Seed_Use <- (Seed_Set[[m]])[SimVal]
  print(paste0("Seed: ",Seed_Use))
  set.seed(Seed_Use)
  Output <- c(Output,
              setNames(list(Analyze_One(Scen=m,
                                 SimNo=SimVal,
                                 NumPerms=Perms_Use,
                                 MVO_list=MVO_list_full,
                                 Comps.Nest=Comp_wts,
                                 CI.GDID=get(paste0("CI.Tx.Obj_",m)),
                                 CI.SV.Mesh=Single_Vals,
                                 CI.Perc=0.95,
                                 Comps=c("TW","CS","MEM","CPI","CPI.T","CPI.D","CPI.DT","CLWP","CLWPA"),
                                 Comps_PermPs=c("TW","CS","MEM","CPI","CPI.T","CPI.D","CPI.DT","CLWP","CLWPA"))),
                       paste0("Sim_",SimVal)))
}
save(Output, file=paste0("sim_res/Scen","_",m,"_ArrNo_",ArrNo,".Rda"))
