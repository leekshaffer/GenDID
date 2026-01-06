#######################################
###### File: Sim_Setup.R ##############
###### Lee Kennedy-Shaffer ############
#######################################

source("R/Analysis.R")
source("R/Sigmas.R")

### Helper Functions for Simulation:

#### Create frame of design:
Sim_Frame <- function(N, J, StartingPds=NULL) {
  if (is.null(StartingPds)) {
    if (J==N+1) {
      warning(simpleWarning(message="Assuming one switch per period starting in period 2"))
      StartingPds <- 2:J
    } else if (J==N) {
      warning(simpleWarning(message="Assuming one switch per period starting in period 1"))
      StartingPds <- 1:J
    } else if (N %% (J-1) == 0) {
      warning(simpleWarning(message=paste0("Assuming ",N/(J-1), " switches per period starting in period 2.")))
      StartingPds <- rep(2:J, each=N/(J-1))
    } else if (N %% J == 0) {
      warning(simpleWarning(message=paste0("Assuming ",N/J, " switches per period starting in period 1.")))
      StartingPds <- rep(1:J, each=N/J)
    } else {
      stop(simpleError(message="N and J must both be numeric. Either N is a multiple of J or J-1, or specify StartingPds."))
    }
  }
  return(tibble(Cluster=rep(1:N, each=J),
                Period=rep(1:J, times=N),
                Start=rep(StartingPds, each=J)) %>%
           dplyr::mutate(Interv=if_else(Period >= Start, 1, 0)))
}

#### Simulate Data:
Sim_Data <- function(Sim.Fr, mu, Alpha1,
                     T1, T2, ProbT1,
                     sig_nu, sig_e, m,
                     ThetaType, ThetaDF) {
  N <- length(unique(Sim.Fr$Cluster))
  J <- length(unique(Sim.Fr$Period))
  Sim.Dat <- Sim.Fr %>%
    dplyr::mutate(FE.g = rep(sample(Alpha1, N, replace=FALSE), each=J),
                  FE.t.Type = rep(sample(c(1,2), size=N, replace=TRUE, prob=c(ProbT1,1-ProbT1)),
                                  each=J),
                  RE.CPI = rnorm(N*J, mean=0, sd=sig_nu),
                  FE.t=if_else(FE.t.Type==1, T1[Period], T2[Period]))
  if (ThetaType==5) {
    Theta.ij <- apply(Sim.Dat, 1,
                      FUN=function(x) x["Interv"]*ThetaDF$Theta[1])
  } else if (ThetaType==4) {
    Theta.ij <- apply(Sim.Dat, 1,
                      FUN=function(x) ifelse(x["Interv"]==0, 0,
                                             ThetaDF[ThetaDF$j==x["Period"],]$Theta))
  } else if (ThetaType==3) {
    Theta.ij <- apply(Sim.Dat, 1,
                      FUN=function(x) ifelse(x["Interv"]==0, 0,
                                             ThetaDF[ThetaDF$a==x["Period"]-x["Start"]+1,]$Theta))
  } else if (ThetaType==2) {
    Theta.ij <- apply(Sim.Dat, 1,
                      FUN=function(x) ifelse(x["Interv"]==0, 0,
                                             ThetaDF[ThetaDF$j==x["Period"] & ThetaDF$a==x["Period"]-x["Start"]+1,]$Theta))
  } else if (ThetaType==1) {
    Theta.ij <- apply(Sim.Dat, 1,
                      FUN=function(x) ifelse(x["Interv"]==0, 0,
                                             ThetaDF[ThetaDF$j==x["Period"] & ThetaDF$a==x["Period"]-x["Start"]+1 & ThetaDF$i==x["Cluster"],]$Theta))
  } else {
    stop(simpleError(message="ThetaType must be an integer from 1 to 5 corresponding to the effect heterogeneity assumption."))
  }
  Sim.Dat <- Sim.Dat %>% dplyr::mutate(Theta.ij=Theta.ij,
                                       mu.ij=mu+FE.g+FE.t+RE.CPI+Theta.ij)
  YVals <- matrix(rnorm(n=N*J*m, mean=rep(Sim.Dat$mu.ij, each=m), sd=sig_e),
                  nrow=N*J, ncol=m, byrow=TRUE)
  colnames(YVals) <- paste("Y.ij", as.character(1:m), sep=".")
  YVals_Out <- as_tibble(YVals) %>%
    mutate(Summ=apply(YVals, 1, FUN=mean),
           Summ_SD=apply(YVals, 1, FUN=sd))
  return(Sim.Dat %>% bind_cols(YVals_Out))
}

### Simulation Scenarios

N_use <- 14
J_use <- 8
mu_use <- 0.3
StartTimes <- tibble(Cluster=1:N_use, StartPd=rep(2:J_use, each=N_use/(J_use-1)))
Sim.Fr <- Sim_Frame(N=N_use, J=J_use, StartingPds=StartTimes$StartPd)

Param_Set <- tibble(
  Scenario=1:9,
  NumSims=rep(1000, 9),
  NumPerms=rep(500, 9),
  mu=rep(mu_use, 9),
  ProbT1=rep(1, 9),
  sig_nu=rep(0.01, 9),
  sig_e=rep(0.1, 9),
  m=100,
  J=J_use,
  N=N_use
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

## Analysis Scenarios

Assns <- 2:5
SigmaNames <- c("CS_0_003","CS_0_333","Ind")

save(list=c("Param_Set",
            "Theta_Set",
            "mu_use",
            "Alpha1",
            "T1",
            "T2",
            "Assns",
            "SigmaNames"),
     file="sim_data/sim-setup.Rda")

## Create and Save Simulation Data

set.seed(73475)

for (i in 1:(dim(Param_Set)[1])) {
  sim_data <- replicate(n=Param_Set$NumSims[i],
                      expr=Sim_Data(Sim.Fr, mu=Param_Set$mu[i], Alpha1=Alpha1,
                                    T1=T1, T2=T2, ProbT1=Param_Set$ProbT1[i],
                                    sig_nu=Param_Set$sig_nu[i], sig_e=Param_Set$sig_e[i], m=Param_Set$m[i],
                                    ThetaType=Theta_Set[[i]]$Type, ThetaDF=Theta_Set[[i]]$ThetaDF),
                      simplify=FALSE)
  assign(x=paste0("sim_data_",i),
         value=sim_data)
  save(list=paste0("sim_data_",i),
       file=paste0("sim_data/sim_data_",i,".Rda"))
}

## Store the Observation Weights for Different Analysis Settings:

### Treatment Effect Vectors:
Vectors.5 <- cbind(Single=c(1))

Vectors.4 <- cbind(AvgEx8=c(rep(1/6,6),0),
                   T.2=c(1,rep(0,6)),
                   T.3=c(0,1,rep(0,5)),
                   T.4=c(0,0,1,rep(0,4)),
                   T.5=c(0,0,0,1,rep(0,3)),
                   T.6=c(rep(0,4),1,0,0),
                   T.7=c(rep(0,5),1,0),
                   T.8=c(rep(0,6),1),
                   Middle=c(0,1/3,1/3,1/3,0,0,0))

Vectors.3 <- cbind(Avg=rep(1/7,7),
                   AvgEx7=c(rep(1/6,6),0),
                   D.1=c(1,rep(0,6)),
                   D.2=c(0,1,rep(0,5)),
                   D.3=c(0,0,1,rep(0,4)),
                   D.4=c(0,0,0,1,rep(0,3)),
                   D.5=c(rep(0,4),1,0,0),
                   D.6=c(rep(0,5),1,0),
                   D.7=c(rep(0,6),1),
                   Middle=c(0,1/3,1/3,1/3,0,0,0))

Vectors.2 <- cbind(Avg=c(rep(1/28,28)),
                   AvgExT8=c(rep(1/21,21),rep(0,7)),
                   D.Avg=1/6*c(1/6,1/6,1/5,1/6,1/5,1/4,1/6,1/5,1/4,1/3,
                               1/6,1/5,1/4,1/3,1/2,1/6,1/5,1/4,1/3,1/2,1,rep(0,7)),
                   T.Avg=1/6*c(1,1/2,1/2,rep(1/3,3),rep(1/4,4),
                               rep(1/5,5),rep(1/6,6),rep(0,7)),
                   D.1=c(1/6,1/6,0,1/6,0,0,1/6,0,0,0,1/6,0,0,0,0,1/6,rep(0,12)),
                   D.2=c(0,0,1/5,0,1/5,0,0,1/5,0,0,0,1/5,0,0,0,0,1/5,rep(0,11)),
                   D.12=c(1/11,1/11,1/11,1/11,1/11,0,1/11,1/11,0,0,1/11,1/11,0,0,0,1/11,1/11,rep(0,11)),
                   T.234=c(rep(1/6,6),rep(0,22)),
                   T.2=c(1,rep(0,27)),
                   T.3=c(rep(0,1),rep(1/2,2),rep(0,25)),
                   T.4=c(rep(0,3),rep(1/3,3),rep(0,22)),
                   T.5=c(rep(0,6),rep(1/4,4),rep(0,18)),
                   T.6=c(rep(0,10),rep(1/5,5),rep(0,13)),
                   T.7=c(rep(0,15),rep(1/6,6),rep(0,7)),
                   U.1=c(1,rep(0,27)),
                   U.2=c(0,1,rep(0,26)),
                   U.3=c(0,0,1,rep(0,25)),
                   U.4=c(rep(0,3),1,rep(0,24)),
                   U.5=c(rep(0,4),1,rep(0,23)),
                   U.6=c(rep(0,5),1,rep(0,22)),
                   U.7=c(rep(0,6),1,rep(0,21)),
                   U.8=c(rep(0,7),1,rep(0,20)),
                   U.9=c(rep(0,8),1,rep(0,19)),
                   U.10=c(rep(0,9),1,rep(0,18)),
                   U.11=c(rep(0,10),1,rep(0,17)),
                   U.12=c(rep(0,11),1,rep(0,16)),
                   U.13=c(rep(0,12),1,rep(0,15)),
                   U.14=c(rep(0,13),1,rep(0,14)),
                   U.15=c(rep(0,14),1,rep(0,13)),
                   U.16=c(rep(0,15),1,rep(0,12)),
                   U.17=c(rep(0,16),1,rep(0,11)),
                   U.18=c(rep(0,17),1,rep(0,10)),
                   U.19=c(rep(0,18),1,rep(0,9)),
                   U.20=c(rep(0,19),1,rep(0,8)),
                   U.21=c(rep(0,20),1,rep(0,7)),
                   Group=1/6*c(1/6,1/5,1/6,1/4,1/5,1/6,1/3,1/4,1/5,1/6,
                               1/2,1/3,1/4,1/5,1/6,1,1/2,1/3,1/4,1/5,1/6,
                               rep(0,7)))

## Variance Assumptions To Use:
SigmaList <- list(Ind=create_Sigma_Ind(N=N_use, J=J_use),
                  CS_0_003=create_Sigma_CS(rho=0.003, N=N_use, J=J_use),
                  CS_0_333=create_Sigma_CS(rho=0.333, N=N_use, J=J_use),
                  AR1_0_012=create_Sigma_AR1(rho=0.012, N=N_use, J=J_use))

MVO_list_full <- NULL
for (k in Assns) {
  SO <- Solve_Assumption(StartTimes,
                         OrderedPds=1:J_use,
                         Assumption=k,
                         v.Mat=get(paste0("Vectors.",k)))

  MVO <- lapply(SigmaList[SigmaNames],
                FUN=function(x) min_var(SolveOut=SO,
                                        Sigma=x,
                                        Drop_ADFT=FALSE,
                                        method="CG", maxit=10000))
  names(MVO) <- paste0("A",k,"_",names(MVO))
  MVO_list_full <- c(MVO_list_full, MVO)
}

## Nested Comparison Approaches Weights

Comp_wts <- Comp_Ests_Weights(ADFT_obj=SO$ADFT,
                              estimator=c("TW","CS","SA","CH","CO","NP"))$Obs.weights

save(list=c("MVO_list_full","Comp_wts"),
     file="sim_data/sim-mvo-list.Rda")


## CI Treatment Effects to Check

Ctrl_prob_overall <- mu_use + mean(Alpha1) + sum(T1*2*(0:7))/sum(2*(0:7))

Single_Vals <- tibble(ValNum=1:701,
                      Name=paste0("SV_",ValNum),
                      Probability=seq(-0.4, 0.3, length.out=701),
                      `Log Odds`=log((Probability+Ctrl_prob_overall)/(1-Probability-Ctrl_prob_overall)) -
                        log(Ctrl_prob_overall/(1-Ctrl_prob_overall)))
SV_list <- split(Single_Vals %>% dplyr::select(Probability,`Log Odds`),
                 seq(nrow(Single_Vals)))
names(SV_list) <- Single_Vals$Name

for (Scen in 1:(dim(Param_Set)[1])) {
  RDs <- Theta_Set[[Scen]]$ThetaDF$Theta
  # CIL <- list(Zero=tibble(Probability=0, `Log Odds`=0))
  if (Theta_Set[[Scen]]$Type == 5) {
    Ctrl_prob <- Ctrl_prob_overall
  } else if (Theta_Set[[Scen]]$Type == 4) {
    Ctrl_prob <- mu_use + mean(Alpha1) + T1[2:8]
  } else if (Theta_Set[[Scen]]$Type == 3) {
    Ctrl_prob <- mu_use + mean(Alpha1) +
      c(mean(T1[2:8]), mean(T1[3:8]), mean(T1[4:8]), mean(T1[5:8]),
        mean(T1[6:8]), mean(T1[7:8]), mean(T1[8:8]))
  }
  CIL <- list(True=tibble(Probability=RDs,
                          `Log Odds`=log((Ctrl_prob+RDs)/(1-Ctrl_prob-RDs)) -
                            log(Ctrl_prob/(1-Ctrl_prob))))
  CIL <- c(CIL, SV_list)

  assign(x=paste0("CI.Tx.Obj_",Scen),
         value=CI_get_Tx_All(ADFT_obj=(MVO_list_full[[paste0("A",Theta_Set[[Scen]]$Type,"_Ind")]])$ADFT,
                             CI_list=CIL))

}

save(list=c("Single_Vals",paste0("CI.Tx.Obj_",(1:(dim(Param_Set)[1])))),
     file="sim_data/sim-CI-objects.Rda")






