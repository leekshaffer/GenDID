#######################################
#### File: Sim_Compile.R ######
#######################################

load("sim_data/sim-setup.Rda")
load("sim_data/sim-CI-objects.Rda")

library(dplyr)
library(tidyr)
library(tibble)

Sim_Results <- NULL
# Missing <- NULL
for (m in Param_Set$Scenario) {
  Scen_Full <- NULL
  for (SimNo in 1:Param_Set$NumSims) {
    if (file.exists(paste0("sim_res/Scen","_",m,"_SimNo_",SimNo,".Rda"))) {
      load(paste0("sim_res/Scen","_",m,"_SimNo_",SimNo,".Rda"))
      Scen_Full <- c(Scen_Full, list(Output))
      rm(Output)
    } else {
      print(paste0("File does not exist for scenario ",m,", simulation number ",SimNo))
      ## To keep list of missing files:
      # Missing <- Missing %>% bind_rows(tibble(Scenario=m, Simulation=SimNo))
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
  Addl_Cvg <- Direct_CI_Vals %>% dplyr::filter(Scenario==m) %>%
    left_join(Frame %>% dplyr::mutate(RowNo=row_number()),
              by=join_by(Estimator,Scenario))
  Addl_Cvg_Vals <- Vals[Addl_Cvg$RowNo,,]
  Addl_Cvg_Mat <- NULL
  for (i in 1:(dim(Addl_Cvg_Vals)[3])) {
    x <- Addl_Cvg_Vals[,,i]
    Addl_Cvg_Mat <- cbind(Addl_Cvg_Mat,x[,"CIL"] <= Addl_Cvg$Estimand & x[,"CIU"] >= Addl_Cvg$Estimand)
  }
  Addl_Cvg <- Addl_Cvg %>% bind_cols(Cvg_New=apply(Addl_Cvg_Mat,1,mean))
  Means[Addl_Cvg$RowNo,"CI Coverage"] <- Addl_Cvg$Cvg_New
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
# save(Missing,
#      file="sim_res/Missing_Simulations.Rda")
save(Sim_Results,
     file="sim_res/Simulation_Results.Rda")

