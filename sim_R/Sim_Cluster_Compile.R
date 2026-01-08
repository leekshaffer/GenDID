#######################################
#### File: Sim_Cluster_Compile.R ######
#######################################

## Make sure this value is set to the same as in Sim_Cluster.R:
NumPerArr <- 10

load("sim_data/sim-setup.Rda")

library(dplyr)
library(tidyr)
library(tibble)

Sim_Results <- NULL
# Errors <- NULL
for (m in Param_Set$Scenario) {
  NumberArrays <- ceiling(Param_Set$NumSims[Param_Set$Scenario==m]/NumPerArr)
  Scen_Full <- NULL
  for (ArrNo in 1:NumberArrays) {
    if (file.exists(paste0("sim_res/Scen","_",m,"_ArrNo_",ArrNo,".Rda"))) {
      load(paste0("sim_res/Scen","_",m,"_ArrNo_",ArrNo,".Rda"))
      Scen_Full <- c(Scen_Full, Output)
      rm(Output)
    } else {
      print(paste0("File does not exist for scenario ",m,", array number ",ArrNo))
      ## To keep list of errors:
      # Errors <- Errors %>% bind_rows(tibble(Scenario=m, Array=ArrNo))
    }
  }
  Frame <- Scen_Full[[1]] %>% dplyr::select(Estimator,Outcome) %>%
    dplyr::mutate(Scenario=m)
  Vals <- sapply(Scen_Full,
                 FUN=function(x) as.matrix(x %>% dplyr::select(-c(Estimator,Outcome))),
                 simplify="array")
  save(list=c("Frame","Vals"),
       file=paste0("sim_res/Sim_Res_Scen_",m,".Rda"))
  Means <- apply(Vals[,c("Estimate","CIW","CI_True","SE"),,drop=FALSE],
                 MARGIN=c(1,2), FUN=mean)
  colnames(Means) <- c("Mean Estimate","Mean CI Width","CI Coverage","Mean SE")
  Medians <- apply(Vals[,"Estimate",,drop=FALSE],
                   MARGIN=c(1,2), FUN=median)
  colnames(Medians) <- "Median Estimate"
  SDs <- apply(Vals[,"Estimate",,drop=FALSE],
               MARGIN=c(1,2), FUN=sd)
  colnames(SDs) <- "SD of Estimate"
  Powers <- apply(Vals[,c("P","P.Asy"),,drop=FALSE],
                  MARGIN=c(1,2), FUN=function(x) mean(x < 0.05))
  colnames(Powers) <- c("Power","Power (Asymptotic)")
  Sim_Results <- Sim_Results %>%
    bind_rows(Frame %>%
                bind_cols(Means,Medians,SDs,Powers))
  ## To check number of simulations:
  print(paste0("Number of Simulations Compiled for Scenario ",m,": ",(dim(Vals)[3])))
  rm(list=c("Scen_Full","Frame","Vals","NumberArrays"))
}
save(Sim_Results,
     file="sim_res/Simulation_Results.Rda")
# save(Errors,
#      file="sim_res/Missing_Jobs.Rda")
