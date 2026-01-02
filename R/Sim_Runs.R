#######################################
###### File: Sim_Runs.R ###############
###### Lee Kennedy-Shaffer ############
#######################################

source("R/Simulations.R")

## A version of simulate_FromSet that allows simple parallelization
simulate_FromSet_Par <- function(Param_Set,
                             Theta_Set,
                             StartingPds=NULL,
                             Alpha1,
                             T1, T2,
                             MVO_list, ADFT_list=NULL,
                             outdir=NULL,
                             outname=NULL,
                             parallel=TRUE, ## allows simple parallelization across settings
                             n_cores=NULL) {
  if (!parallel) {
    simulate_FromSet(Param_Set, Theta_Set, StartingPds,
                     Alpha1, T1, T2,
                     MVO_list, ADFT_list,
                     outdir, outname)
  } else {
    require(foreach)
    require(doParallel)
    if (is.null(n_cores)) {
      n_cores <- detectCores() - 1
    }
    par_clust <- makeCluster(n_cores)
    registerDoParallel(par_clust)
    foreach (i=1:(dim(Param_Set)[1])) %dopar% {
      source("R/Simulations.R")
      row <- Param_Set[i,]
      print(paste("Starting Sim. Setting Number",row$SimNo))
      assign(x=paste0("Res_Sim_",i),
             value=simulate_SWT(row$NumSims,
                                row$N, row$J, StartingPds,
                                row$mu, Alpha1,
                                T1, T2, row$ProbT1,
                                row$sig_nu, row$sig_e, row$m,
                                Theta_Set[[i]]$Type, Theta_Set[[i]]$ThetaDF,
                                MVO_list, ADFT_list,
                                Comparisons=Theta_Set[[i]]$Comps, Theta_Set[[i]]$corstr,
                                Permutations=row$NumPerms))
      save(list=paste0("Res_Sim_",row$SimNo),
           file=paste0(outdir,"/",outname,"_",row$SimNo,".Rda"))
    }
    stopCluster(cl = par_clust)
  }
}

## Simulation runs:

### Parameters for Simulation:
NumSims.all <- 10 # Full: 1000
NumPerms.all <- 50 # Full: 250
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

Theta_Set <- list(list(Type=5, ThetaDF=tibble(Theta=0), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.D","CPI.DT","CLWP","CLWPA"), corstr=NULL),
                  list(Type=5, ThetaDF=tibble(Theta=-0.02), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.D","CPI.DT","CLWP","CLWPA"), corstr=NULL),
                  list(Type=5, ThetaDF=tibble(Theta=-0.04), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.D","CPI.DT","CLWP","CLWPA"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=seq(from=-0.07,to=0.05,by=0.02)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.DT","CLWP","CLWPA"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=c(-0.07,-0.06,-0.04,0,0.03,0.02,0.01)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.DT","CLWP","CLWPA"), corstr=NULL),
                  list(Type=4, ThetaDF=tibble(j=2:8, Theta=c(rep(-0.03,4),rep(0,3))), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.T","CPI.DT","CLWP","CLWPA"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=seq(from=-0.01,to=-0.04,by=-0.005)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D","CPI.DT","CLWP","CLWPA"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=c(rep(0,2),rep(-0.03,5))), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D","CPI.DT","CLWP","CLWPA"), corstr=NULL),
                  list(Type=3, ThetaDF=tibble(a=1:7, Theta=seq(from=-0.07,to=0.05,by=0.02)), Comps=c("TW","CS","SA","CH","CO","NP","CPI","CPI.D","CPI.DT","CLWP","CLWPA"), corstr=NULL))

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

## Timing:
### 3 cores on home comp, 10 sims/scenario, 25 perms/sim: 15 mins
### 3 cores on home comp, 10 sims/scenario, 50 perms/sim: 15 mins
pt <- proc.time()

set.seed(73475)

folder <- "sim_res_new" ## sim_res

simulate_FromSet_Par(Param_Set, Theta_Set,
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
                 # ADFT_list=list(Comp=ADFT_5),
                 ADFT_list=NULL,
                 outdir=folder, outname="Sim_Set",
                 parallel=TRUE,
                 n_cores=9)



### Get simulation results:
Full_Sim_Res <- NULL
for (i in Param_Set$SimNo) {
  load(paste0(folder,"/Sim_Set_",i,".Rda"))
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
     file=paste0(folder,"/Full_Sim_Res.Rda"))

proc.time() - pt
