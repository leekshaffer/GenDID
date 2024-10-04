#######################################
###### File: Sim_Runs.R ###############
###### Lee Kennedy-Shaffer ############
###### Created 2024/10/04 #############
#######################################

source("Simulations.R")

## Original simulation runs (1--2):

### Parameters for Simulation:
NumSims.all <- 250
NumPerms.all <- 250
Param_Set <- tribble(
  ~SimNo, ~NumSims, ~NumPerms, ~mu, ~ProbT1, ~sig_nu, ~sig_e, ~m, ~J, ~N,
  1, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  2, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  3, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  10, NumSims.all, NumPerms.all, 0.3, 0.5, 0.01, 0.1, 100, 8, 14,
  11, NumSims.all, NumPerms.all, 0.3, 0.5, 0.01, 0.1, 100, 8, 14,
  12, NumSims.all, NumPerms.all, 0.3, 0.5, 0.01, 0.1, 100, 8, 14,
  4, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  5, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  6, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  7, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  8, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  9, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  13, NumSims.all, NumPerms.all, 0.3, 0.5, 0.01, 0.1, 100, 8, 14,
  14, NumSims.all, NumPerms.all, 0.3, 0.5, 0.01, 0.1, 100, 8, 14,
  15, NumSims.all, NumPerms.all, 0.3, 0.5, 0.01, 0.1, 100, 8, 14,
  16, NumSims.all, NumPerms.all, 0.3, 0.5, 0.01, 0.1, 100, 8, 14,
  17, NumSims.all, NumPerms.all, 0.3, 0.5, 0.01, 0.1, 100, 8, 14,
  18, NumSims.all, NumPerms.all, 0.3, 0.5, 0.01, 0.1, 100, 8, 14
)

Theta_Set <- list(list(Type=5, ThetaDF=tibble(Theta=0), Comps=c("TW","CS","SA","CH","CO","NP","CPI"), corstr=NULL),
                  list(Type=5, ThetaDF=tibble(Theta=-0.02), Comps=c("TW","CS","SA","CH","CO","NP","CPI"), corstr=NULL),
                  list(Type=5, ThetaDF=tibble(Theta=-0.04), Comps=c("TW","CS","SA","CH","CO","NP","CPI"), corstr=NULL),
                  list(Type=5, ThetaDF=tibble(Theta=0), Comps=c("TW","CS","SA","CH","CO","NP","CPI"), corstr=NULL),
                  list(Type=5, ThetaDF=tibble(Theta=-0.02), Comps=c("TW","CS","SA","CH","CO","NP","CPI"), corstr=NULL),
                  list(Type=5, ThetaDF=tibble(Theta=-0.04), Comps=c("TW","CS","SA","CH","CO","NP","CPI"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=seq(from=-0.07,to=0.05,by=0.02)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=c(-0.07,-0.06,-0.04,0,0.03,0.02,0.01)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=c(rep(-0.03,4),rep(0,3))), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=seq(from=-0.01,to=-0.04,by=-0.005)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=c(rep(0,2),rep(-0.03,5))), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=seq(from=-0.07,to=0.05,by=0.02)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=seq(from=-0.07,to=0.05,by=0.02)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=c(-0.07,-0.06,-0.04,0,0.03,0.02,0.01)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=c(rep(-0.03,4),rep(0,3))), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=seq(from=-0.01,to=-0.04,by=-0.005)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=c(rep(0,2),rep(-0.03,5))), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=seq(from=-0.07,to=0.05,by=0.02)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D"), corstr=NULL))

Alpha1 <- c(-0.007, 0.003, 0.008, -0.016, -0.003, -0.005, -0.012, 
            0.002, 0.005, -0.001, 0.020, 0, 0.017, -0.011)
T1 <- c(0,0.08,0.18,0.29,0.30,0.27,0.20,0.13)
T2 <- c(0,0.02,0.03,0.07,0.13,0.19,0.27,0.30)

### Run Simulations:

load("int/xpert-mv-a_2_CS_0_003.Rda")
load("int/xpert-mv-a_3_CS_0_003.Rda")
load("int/xpert-mv-a_4_CS_0_003.Rda")
load("int/xpert-mv-a_5_CS_0_003.Rda")
load("int/xpert-mv-a_2_CS_0_333.Rda")
load("int/xpert-mv-a_3_CS_0_333.Rda")
load("int/xpert-mv-a_4_CS_0_333.Rda")
load("int/xpert-mv-a_5_CS_0_333.Rda")
load("int/xpert-mv-a_2_Ind.Rda")
load("int/xpert-mv-a_3_Ind.Rda")
load("int/xpert-mv-a_4_Ind.Rda")
load("int/xpert-mv-a_5_Ind.Rda")
load("../int_large/xpert-solve-a_2.Rda")
load("../int_large/xpert-solve-a_3.Rda")
load("../int_large/xpert-solve-a_4.Rda")
load("../int_large/xpert-solve-a_5.Rda")

## Run 1 seed:
set.seed(801611)
## Run 2 seed:
set.seed(569908)

system.time(simulate_FromSet(Param_Set,
                             Theta_Set,
                             StartingPds=NULL,
                             Alpha1, 
                             T1, T2, 
                             MVO_list=list(A2_003=MVOut_2_CS_0_003,
                                           A3_003=MVOut_3_CS_0_003,
                                           A4_003=MVOut_4_CS_0_003,
                                           A5_003=MVOut_5_CS_0_003,
                                           A2_333=MVOut_2_CS_0_333,
                                           A3_333=MVOut_3_CS_0_333,
                                           A4_333=MVOut_4_CS_0_333,
                                           A5_333=MVOut_5_CS_0_333,
                                           A2_Ind=MVOut_2_Ind,
                                           A3_Ind=MVOut_3_Ind,
                                           A4_Ind=MVOut_4_Ind,
                                           A5_Ind=MVOut_5_Ind),
                             SO_list=list(Comp=SolveOut_5),
                             outdir="sim_res_2", # change to sim_res_1 for run 1
                             outname="Sim_Set"))



### Get simulation results:
Full_Sim_Res <- NULL
for (i in Param_Set$SimNo) {
  load(paste0("sim_res_1/Sim_Set_",i,".Rda"))
  assign(x="Res1", value=get(paste0("Res_Sim_",i)))
  load(paste0("sim_res_2/Sim_Set_",i,".Rda"))
  assign(x="Res2", value=get(paste0("Res_Sim_",i)))
  Ests <- rbind(Res1$Estimates, Res2$Estimates)
  PVs <- rbind(Res1$PValues, Res2$PValues)
  Res_int <- tibble(SimNo=rep(i,4), 
                    Result=c("Mean Estimate","Median Estimate","SD Estimate","Power"))
  Full_Sim_Res <- Full_Sim_Res %>% bind_rows(cbind(Res_int, rbind(apply(Ests, 2, mean, na.rm=TRUE),
                                                                  apply(Ests, 2, median, na.rm=TRUE),
                                                                  apply(Ests, 2, sd, na.rm=TRUE),
                                                                  apply(PVs, 2, function(x) mean(x <= 0.05, na.rm=TRUE)))))
  rm(list=c("Res1","Res2","Ests","PVs",
            paste0("Res_Sim_",i),"Res_int"))
}
save(Full_Sim_Res,
     file="sim_res/Full_Sim_Res_12.Rda")
write_csv(x=Full_Sim_Res,
          file="sim_res/Full_Sim_Res_12.csv")



## New simulation runs (3):

### Parameters for Simulation:
NumSims.all <- 500
NumPerms.all <- 250
Param_Set <- tribble(
  ~SimNo, ~NumSims, ~NumPerms, ~mu, ~ProbT1, ~sig_nu, ~sig_e, ~m, ~J, ~N,
  1, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  2, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  3, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  4, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  5, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  6, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  7, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  8, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14,
  9, NumSims.all, NumPerms.all, 0.3, 1, 0.01, 0.1, 100, 8, 14
)

Theta_Set <- list(list(Type=5, ThetaDF=tibble(Theta=0), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.D","CPI.DT"), corstr=NULL),
                  list(Type=5, ThetaDF=tibble(Theta=-0.02), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.D","CPI.DT"), corstr=NULL),
                  list(Type=5, ThetaDF=tibble(Theta=-0.04), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.D","CPI.DT"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=seq(from=-0.07,to=0.05,by=0.02)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.DT"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=c(-0.07,-0.06,-0.04,0,0.03,0.02,0.01)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.DT"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=c(rep(-0.03,4),rep(0,3))), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.DT"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=seq(from=-0.01,to=-0.04,by=-0.005)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D","CPI.DT"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=c(rep(0,2),rep(-0.03,5))), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D","CPI.DT"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=seq(from=-0.07,to=0.05,by=0.02)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D","CPI.DT"), corstr=NULL))

Alpha1 <- c(-0.007, 0.003, 0.008, -0.016, -0.003, -0.005, -0.012, 
            0.002, 0.005, -0.001, 0.020, 0, 0.017, -0.011)
T1 <- c(0,0.08,0.18,0.29,0.30,0.27,0.20,0.13)
T2 <- c(0,0.02,0.03,0.07,0.13,0.19,0.27,0.30)

### Run Simulations:

load("int/xpert-mv-a_2_CS_0_003.Rda")
load("int/xpert-mv-a_3_CS_0_003.Rda")
load("int/xpert-mv-a_4_CS_0_003.Rda")
load("int/xpert-mv-a_5_CS_0_003.Rda")
load("int/xpert-mv-a_2_CS_0_333.Rda")
load("int/xpert-mv-a_3_CS_0_333.Rda")
load("int/xpert-mv-a_4_CS_0_333.Rda")
load("int/xpert-mv-a_5_CS_0_333.Rda")
load("int/xpert-mv-a_2_Ind.Rda")
load("int/xpert-mv-a_3_Ind.Rda")
load("int/xpert-mv-a_4_Ind.Rda")
load("int/xpert-mv-a_5_Ind.Rda")
load("../int_large/xpert-solve-a_2.Rda")
load("../int_large/xpert-solve-a_3.Rda")
load("../int_large/xpert-solve-a_4.Rda")
load("../int_large/xpert-solve-a_5.Rda")

set.seed(73475)

simulate_FromSet(Param_Set, Theta_Set,
                 StartingPds=NULL, 
                 Alpha1,T1, T2, 
                 MVO_list=list(A2_003=MVOut_2_CS_0_003,
                               A3_003=MVOut_3_CS_0_003,
                               A4_003=MVOut_4_CS_0_003,
                               A5_003=MVOut_5_CS_0_003,
                               A2_333=MVOut_2_CS_0_333,
                               A3_333=MVOut_3_CS_0_333,
                               A4_333=MVOut_4_CS_0_333,
                               A5_333=MVOut_5_CS_0_333,
                               A2_Ind=MVOut_2_Ind,
                               A3_Ind=MVOut_3_Ind,
                               A4_Ind=MVOut_4_Ind,
                               A5_Ind=MVOut_5_Ind),
                 SO_list=list(Comp=SolveOut_5),
                 outdir="sim_res_3", outname="Sim_Set")



### Get simulation results:
Full_Sim_Res <- NULL
for (i in Param_Set$SimNo) {
  load(paste0("sim_res_3/Sim_Set_",i,".Rda"))
  assign(x="Res1", value=get(paste0("Res_Sim_",i)))
  Ests <- Res1$Estimates
  PVs <- Res1$PValues
  Res_int <- tibble(SimNo=rep(i,4),
                    Result=c("Mean Estimate","Median Estimate","SD Estimate","Power"))
  Full_Sim_Res <- Full_Sim_Res %>% bind_rows(cbind(Res_int, rbind(apply(Ests, 2, mean, na.rm=TRUE),
                                                                  apply(Ests, 2, median, na.rm=TRUE),
                                                                  apply(Ests, 2, sd, na.rm=TRUE),
                                                                  apply(PVs, 2, function(x) mean(x <= 0.05, na.rm=TRUE)))))
  rm(list=c("Res1","Ests","PVs",
            paste0("Res_Sim_",i),"Res_int"))
}
save(Full_Sim_Res,
     file="sim_res/Full_Sim_Res_3.Rda")
write_csv(x=Full_Sim_Res,
          file="sim_res/Full_Sim_Res_3.csv")
