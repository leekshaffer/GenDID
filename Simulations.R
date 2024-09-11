#######################################
###### File: Simulations.R ###########
###### Lee Kennedy-Shaffer ############
###### Created 2024/09/09 #############
#######################################

source("A_Const.R")
source("Sigmas.R")
source("Full_Analysis.R")
source("CompEsts.R")

require(tidyverse)
require(lme4)
require(geepack)

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
      warning(simpleWarning(message=paste0("Assuming ",N/(J-1), " switches per period starting in period 1.")))
      StartingPds <- rep(2:J, each=N/(J-1))
    } else if (N %% J == 0) {
      warning(simpleWarning(message=paste0("Assuming ",N/J, " switches per period starting in period 1.")))
      StartingPds <- rep(1:J, each=N/J)
    } else {
      stop(simpleError(message="N and J must both be numeric."))
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
                     sig_nu, sig_e, Theta, m) {
  Sim.Dat <- Sim.Fr %>%
    dplyr::mutate(FE.g = rep(sample(Alpha1, N, replace=FALSE), each=J),
                  FE.t.Type = rep(sample(c(1,2), size=N, replace=TRUE, prob=c(ProbT1,1-ProbT1)), 
                                  each=J),
                  CPI = rnorm(N*J, mean=0, sd=sig_nu),
                  FE.t=if_else(FE.t.Type==1, T1[Period], T2[Period]),
                  mu.ij=mu+FE.g+FE.t+CPI+Interv*Theta)
  YVals <- matrix(rnorm(n=N*J*m, mean=rep(Sim.Dat$mu.ij, each=m), sd=sig_e),
                  nrow=N*J, ncol=m, byrow=TRUE)
  colnames(YVals) <- paste("Y.ij", as.character(1:m), sep=".")
  YVals_Out <- as_tibble(YVals) %>%
    mutate(Y.ij.bar=apply(YVals, 1, FUN=mean),
           Y.ij.sd=apply(YVals, 1, FUN=sd))
  return(Sim.Dat %>% bind_cols(YVals_Out))
}

### Analyze Simulated Data
#### (run after Solve_Assumption & MV_Assumption (w/ no Obs or Perms) & Comp_Ests_Weights for setting):
#### Note MV_Out can be a list with multiple MV_Out objects: names are passed on
#### Solve_Out is used only for comparisons to other estimators; can be a list; names passed on
Sim_Weights <- function(MV_Out,
                        Solve_Out=NULL, 
                        Comparisons=NULL) {
  Weights <- NULL
  if (!is.null(MV_Out$Obs.Weights)) { #If only one MV_Out object is given
    Weights <- tibble(GenDID=MV_Out$Obs.Weights)
  } else {
    for (i in 1:length(MV_Out)) { #If list of MV_Out objects is given
      Weights <- Weights %>% 
        bind_cols(as_tibble(MV_Out[[i]]$MV$Obs.weights) %>% 
                    dplyr::rename_with(~paste(names(MV_Out)[i],colnames(MV_Out[[i]]$MV$Obs.weights),
                                              sep="_",recycle0=TRUE)))
    }
  }
  CE_Weights <- NULL
  if (sum(c("TW","CS","SA","CH","CO","NP") %in% Comparisons) > 0) {
    if (!is.null(Solve_Out$DFT)) { #If only one Solve_Out object is given
      CE_Weights <- Comp_Ests_Weights(DFT_obj=Solve_Out$DFT,
                                      Amat=Solve_Out$Amat,
                                      estimator=Comparisons[Comparisons %in% c("TW","CS","SA","CH","CO","NP")])$Obs.weights
    } else {
      for (i in 1:length(Solve_Out)) { #If list of Solve_Out objects is given
        CEwi <- Comp_Ests_Weights(DFT_obj=Solve_Out[[i]]$DFT,
                                  Amat=Solve_Out[[i]]$Amat,
                                  estimator=Comparisons[Comparisons %in% c("TW","CS","SA","CH","CO","NP")])$Obs.weights
        colnames(CEwi) <- paste(names(Solve_Out)[i], colnames(CEwi), sep="_", recycle0=TRUE)
        CE_Weights <- CE_Weights %>% 
          bind_cols(CEwi)
      }
    }
  }
  return(Weights %>% bind_cols(CE_Weights))
}

Sim_Analyze <- function(Sim.Dat,
                        Sim.Wt,
                        MEM=FALSE,
                        CPI=FALSE,
                        GEE=FALSE, #Note: GEE comp is very slow; ~1min/GEE fit
                        corstr="exchangeable") {
  Obs.mat<- t(Sim.Dat[,"Y.ij.bar",,drop=TRUE])
  Results <- Obs.mat %*% as.matrix(Sim.Wt)
  if (MEM | CPI | GEE) {
    Sim.Dat.long <- apply(Sim.Dat, 3, 
                          FUN=function(x) as_tibble(x) %>% 
                            dplyr::select(-c("Y.ij.bar","Y.ij.sd")) %>%
                            pivot_longer(cols=starts_with("Y.ij."),
                                         names_to="Indiv",
                                         names_prefix="Y.ij.",
                                         values_to="Y"))
    if (MEM) {
      MEM_Res <- sapply(Sim.Dat.long,
                    FUN=function(x) unname(fixef(lmer(Y~Interv+factor(Period)+(1|Cluster), data=x))["Interv"]))
      Results <- cbind(Results, MEM=MEM_Res)
    }
    if (CPI) {
      CPI_Res <- sapply(Sim.Dat.long,
                    FUN=function(x) unname(fixef(lmer(Y~Interv+factor(Period)+(1|Cluster)+(1|CPI), 
                                                      data=x %>% dplyr::mutate(CPI=paste(Cluster,Period,sep="_"))))["Interv"]))
      Results <- cbind(Results, CPI=CPI_Res)
    }
    if (GEE) {
      GEE_Res <- sapply(Sim.Dat.long,
                    FUN=function(x) unname(geeglm(Y~Interv+factor(Period),
                                                  data=x %>% dplyr::mutate(Clf=factor(Cluster)),
                                                  id=Clf,
                                                  corstr=corstr)$coefficients["Interv"]))
      Results <- cbind(Results, GEE=GEE_Res)
    }
  }
  return(Results)
}

Sim_Permutations <- function(Sim.Dat,
                        Sim.Wt,
                        Permutations,
                        N=NULL, J=NULL,
                        MEM=FALSE,
                        CPI=FALSE,
                        GEE=FALSE, #Note: GEE comp is very slow
                        corstr="exchangeable") {
  if (is.null(N) | is.null(J)) {
    warning(simpleWarning(message="Getting N and J from Sim.Dat"))
    N <- length(unique(Sim.Dat[,"Cluster",1,drop=TRUE]))
    J <- length(unique(Sim.Dat[,"Period",1,drop=TRUE]))
  }
  Orders <- replicate(n=dim(Sim.Dat)[3],
                      expr=Permute_order(N, J))
  Obs.true <- Sim.Dat[,"Y.ij.bar",,drop=TRUE]
  Perm.mat <- t(sapply(1:dim(Sim.Dat)[3],
                     FUN=function(i) Obs.true[Orders[,i],i]))
  Res.perm <- Perm.mat %*% as.matrix(Sim.Wt)
}

### Parameters for Simulation:
mu <- 0.3
Alpha1 <- c(-0.007, 0.003, 0.008, -0.016, -0.003, -0.005, -0.012, 
            0.002, 0.005, -0.001, 0.020, 0, 0.017, -0.011)
T1 <- c(0,0.08,0.18,0.29,0.30,0.27,0.20,0.13)
T2 <- c(0,0.02,0.03,0.07,0.13,0.19,0.27,0.30)
ProbT1 <- 1
sig_nu <- 0.01
sig_e <- 0.1
m <- 100
J <- 8
N <- 14
Theta <- 0

### Framework for simulated data with fixed treatment schedule:
simulate_SWT <- function(NumSims,
                          N, J, StartingPds=NULL,
                          mu, Alpha1, 
                          T1, T2, ProbT1,
                          sig_nu, sig_e, Theta, m) {
  Sim.Fr <- Sim_Frame(N, J, StartingPds)
  return(replicate(n=NumSims, as.matrix(Sim_Data(Sim.Fr, mu, Alpha1, 
                                                 T1, T2, ProbT1,
                                                 sig_nu, sig_e, Theta, m)), 
                   simplify="array"))
}

### Generate data:
set.seed(801611)
system.time(a <- simulate_SWT(100, N, J, StartingPds=NULL,
                              mu, Alpha1, T1, T2, ProbT1, sig_nu, sig_e, Theta, m))


