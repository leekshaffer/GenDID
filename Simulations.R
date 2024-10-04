#######################################
###### File: Simulations.R ###########
###### Lee Kennedy-Shaffer ############
###### Created 2024/09/09 #############
#######################################

## Note: some details of the mixed effects model fitting is specific to the setting used in Sim_Runs.R
### Please check and alter these or do not use those models if using a different setup.

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
                                              sep="_",recycle0=FALSE)))
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
        colnames(CEwi) <- paste(names(Solve_Out)[i], colnames(CEwi), sep="_", recycle0=FALSE)
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
                        CPI.T=FALSE,
                        CPI.D=FALSE,
                        CPI.DT=FALSE,
                        GEE=FALSE, #Note: GEE comp is very slow; ~1min/GEE fit
                        corstr="exchangeable") {
  Obs.mat<- t(Sim.Dat[,"Y.ij.bar",,drop=TRUE])
  Results <- Obs.mat %*% as.matrix(Sim.Wt)
  if (MEM | CPI | CPI.T | CPI.D | CPI.DT | GEE) {
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
      Results <- cbind(Results, Comp_MEM=MEM_Res)
    }
    if (CPI) {
      CPI_Res <- sapply(Sim.Dat.long,
                    FUN=function(x) unname(fixef(lmer(Y~Interv+factor(Period)+(1|Cluster)+(1|CPI), 
                                                      data=x %>% dplyr::mutate(CPI=paste(Cluster,Period,sep="_"))))["Interv"]))
      Results <- cbind(Results, Comp_CPI=CPI_Res)
    }
    if (CPI.T) {
      CPI.T_Res <- as_tibble(t(sapply(Sim.Dat.long,
                          FUN=function(x) fixef(lmer(Y~Interv*factor(Period) + (1|Cluster) + (1|CPI),
                                                            data=x %>% dplyr::mutate(CPI=paste(Cluster,Period,sep="_")),
                                                     control=lmerControl(check.rankX="silent.drop.cols",
                                                                         check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))))) %>%
        dplyr::select(starts_with("Interv")) %>%
        dplyr::mutate(across(.cols=-c("Interv"),
                             .fns=~.x+Interv)) %>%
        dplyr::rename(`Interv:factor(Period)7`=Interv)
      colnames(CPI.T_Res) <- gsub(".*)","CPI.T_",colnames(CPI.T_Res))
      CPI.T_Res$CPI.T_AvgExLast <- apply(CPI.T_Res, 1, mean, na.rm=TRUE)
      CPI.T_Res$CPI.T_Middle <- apply(CPI.T_Res[,paste0("CPI.T_",2:4)], 1, mean, na.rm=TRUE)
      Results <- cbind(Results, as.matrix(CPI.T_Res))
    }
    if (CPI.D) {
      CPI.D_Res <- as_tibble(t(sapply(Sim.Dat.long,
                                      FUN=function(x) fixef(lmer(Y~Interv+Interv:factor(Diff) + factor(Period)+ (1|Cluster) + (1|CPI),
                                                                 data=x %>% dplyr::mutate(CPI=paste(Cluster,Period,sep="_"),
                                                                                          Diff=if_else(Period-Start < 0,0,Period-Start+1)),
                                                                 control=lmerControl(check.rankX="silent.drop.cols",
                                                                                     check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))))) %>%
        dplyr::select(starts_with("Interv")) %>%
        dplyr::mutate(across(.cols=-c("Interv"),
                             .fns=~.x+Interv)) %>%
        dplyr::rename(`Interv:factor(Diff)7`=Interv)
      colnames(CPI.D_Res) <- gsub(".*)","CPI.D_",colnames(CPI.D_Res))
      CPI.D_Res$CPI.D_Avg <- apply(CPI.D_Res, 1, mean, na.rm=TRUE)
      CPI.D_Res$CPI.D_Middle <- apply(CPI.D_Res[,paste0("CPI.D_",2:4)], 1, mean, na.rm=TRUE)
      Results <- cbind(Results, as.matrix(CPI.D_Res))
    }
    if (CPI.DT) {
      CPI.DT_Res <- as_tibble(t(sapply(Sim.Dat.long,
                                      FUN=function(x) fixef(lmer(Y~Interv+Interv:factor(Diff):factor(Period) + 
                                                                   factor(Period)+ (1|Cluster) + (1|CPI),
                                                                 data=x %>% dplyr::mutate(CPI=paste(Cluster,Period,sep="_"),
                                                                                          Diff=if_else(Period-Start < 0,0,Period-Start+1)),
                                                                 control=lmerControl(check.rankX="silent.drop.cols",
                                                                                     check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))))) %>%
        dplyr::select(starts_with("Interv")) %>%
        dplyr::mutate(across(.cols=-c("Interv"),
                             .fns=~.x+Interv)) %>%
        dplyr::rename(`Interv:factor(Diff)6:factor(Period)7`=Interv)
      colnames(CPI.DT_Res) <- paste0("CPI.DT_a",substr(colnames(CPI.DT_Res),20,20),
                                     "_j",substr(colnames(CPI.DT_Res),36,36))
      CPI.DT_Res_orig <- CPI.DT_Res
      CPI.DT_Res$CPI.DT_Avg <- apply(CPI.DT_Res_orig, 1, mean, na.rm=TRUE)
      CPI.DT_Res$CPI.DT_AvgEx8 <- apply(CPI.DT_Res_orig %>% dplyr::select(-ends_with("j8")), 1, mean, na.rm=TRUE)
      for (i in 1:6) {
        CPI.DT_Res[,paste0("CPI.DT_D",i)] <- apply(CPI.DT_Res_orig %>% dplyr::select(contains(paste0("a",i))), 
                                                   1, mean, na.rm=TRUE)
      }
      for (i in 2:7) {
        CPI.DT_Res[,paste0("CPI.DT_T",i)] <- apply(CPI.DT_Res_orig %>% dplyr::select(contains(paste0("j",i))), 
                                                   1, mean, na.rm=TRUE)
      }
      CPI.DT_Res$CPI.DT_DAvg <- apply(CPI.DT_Res %>% dplyr::select(starts_with("CPI.DT_D")),
                                      1, mean, na.rm=TRUE)
      CPI.DT_Res$CPI.DT_TAvg <- apply(CPI.DT_Res %>% dplyr::select(starts_with("CPI.DT_T")),
                                      1, mean, na.rm=TRUE)
      Results <- cbind(Results, as.matrix(CPI.DT_Res))
    }
    if (GEE) {
      GEE_Res <- sapply(Sim.Dat.long,
                    FUN=function(x) unname(geeglm(Y~Interv+factor(Period),
                                                  data=x %>% dplyr::mutate(Clf=factor(Cluster)),
                                                  id=Clf,
                                                  corstr=corstr)$coefficients["Interv"]))
      Results <- cbind(Results, Comp_GEE=GEE_Res)
    }
  }
  return(Results)
}

Permute_All <- function(Sim.Dat.1,
                        Order1) {
  return(as_tibble(Sim.Dat.1) %>% 
           dplyr::select(all_of(c("Cluster","Period","Start","Interv"))) %>%
           bind_cols((as_tibble(Sim.Dat.1) %>% 
                        dplyr::select(starts_with("Y.ij.")))[Order1,]) %>%
           dplyr::select(-c("Y.ij.bar","Y.ij.sd")) %>%
           pivot_longer(cols=starts_with("Y.ij."),
                        names_to="Indiv",
                        names_prefix="Y.ij.",
                        values_to="Y"))
}

Sim_Permutation <- function(Sim.Dat,
                        Sim.Wt,
                        N=NULL, J=NULL,
                        MEM=FALSE,
                        CPI=FALSE,
                        CPI.T=FALSE,
                        CPI.D=FALSE,
                        CPI.DT=FALSE,
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
  if (MEM | CPI | CPI.T | CPI.D | CPI.DT | GEE) {
    Perm.Dat.Long <- lapply(1:(dim(Sim.Dat)[3]),
                            FUN=function(i) Permute_All(Sim.Dat[,,i], Orders[,i]))
    if (MEM) {
      MEM_Res <- sapply(Perm.Dat.Long,
                        FUN=function(x) unname(fixef(lmer(Y~Interv+factor(Period)+(1|Cluster), data=x))["Interv"]))
      Res.perm <- cbind(Res.perm, Comp_MEM=MEM_Res)
    }
    if (CPI) {
      CPI_Res <- sapply(Perm.Dat.Long,
                        FUN=function(x) unname(fixef(lmer(Y~Interv+factor(Period)+(1|Cluster)+(1|CPI), 
                                                          data=x %>% dplyr::mutate(CPI=paste(Cluster,Period,sep="_"))))["Interv"]))
      Res.perm <- cbind(Res.perm, Comp_CPI=CPI_Res)
    }
    if (CPI.T) {
      CPI.T_Res <- as_tibble(t(sapply(Perm.Dat.Long,
                                      FUN=function(x) fixef(lmer(Y~Interv*factor(Period) + (1|Cluster) + (1|CPI),
                                                                 data=x %>% dplyr::mutate(CPI=paste(Cluster,Period,sep="_")),
                                                                 control=lmerControl(check.rankX="silent.drop.cols",
                                                                                     check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))))) %>%
        dplyr::select(starts_with("Interv")) %>%
        dplyr::mutate(across(.cols=-c("Interv"),
                             .fns=~.x+Interv)) %>%
        dplyr::rename(`Interv:factor(Period)0`=Interv)
      colnames(CPI.T_Res) <- gsub(".*)","CPI.T_",colnames(CPI.T_Res))
      CPI.T_Res$CPI.T_AvgExLast <- apply(CPI.T_Res, 1, mean, na.rm=TRUE)
      CPI.T_Res$CPI.T_Middle <- apply(CPI.T_Res[,paste0("CPI.T_",2:4)], 1, mean, na.rm=TRUE)
      Res.perm <- cbind(Res.perm, as.matrix(CPI.T_Res))
    }
    if (CPI.D) {
      CPI.D_Res <- as_tibble(t(sapply(Perm.Dat.Long,
                                      FUN=function(x) fixef(lmer(Y~Interv+Interv:factor(Diff) + factor(Period)+ (1|Cluster) + (1|CPI),
                                                                 data=x %>% dplyr::mutate(CPI=paste(Cluster,Period,sep="_"),
                                                                                          Diff=if_else(Period-Start < 0,0,Period-Start+1)),
                                                                 control=lmerControl(check.rankX="silent.drop.cols",
                                                                                     check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))))) %>%
        dplyr::select(starts_with("Interv")) %>%
        dplyr::mutate(across(.cols=-c("Interv"),
                             .fns=~.x+Interv)) %>%
        dplyr::rename(`Interv:factor(Diff)0`=Interv)
      colnames(CPI.D_Res) <- gsub(".*)","CPI.D_",colnames(CPI.D_Res))
      CPI.D_Res$CPI.D_Avg <- apply(CPI.D_Res, 1, mean, na.rm=TRUE)
      CPI.D_Res$CPI.D_Middle <- apply(CPI.D_Res[,paste0("CPI.D_",2:4)], 1, mean, na.rm=TRUE)
      Res.perm <- cbind(Res.perm, as.matrix(CPI.D_Res))
    }
    if (CPI.DT) {
      CPI.DT_Res <- as_tibble(t(sapply(Perm.Dat.Long,
                                       FUN=function(x) fixef(lmer(Y~Interv+Interv:factor(Diff):factor(Period) + 
                                                                    factor(Period)+ (1|Cluster) + (1|CPI),
                                                                  data=x %>% dplyr::mutate(CPI=paste(Cluster,Period,sep="_"),
                                                                                           Diff=if_else(Period-Start < 0,0,Period-Start+1)),
                                                                  control=lmerControl(check.rankX="silent.drop.cols",
                                                                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))))) %>%
        dplyr::select(starts_with("Interv")) %>%
        dplyr::mutate(across(.cols=-c("Interv"),
                             .fns=~.x+Interv)) %>%
        dplyr::rename(`Interv:factor(Diff)6:factor(Period)7`=Interv)
      colnames(CPI.DT_Res) <- paste0("CPI.DT_a",substr(colnames(CPI.DT_Res),20,20),
                                     "_j",substr(colnames(CPI.DT_Res),36,36))
      CPI.DT_Res_orig <- CPI.DT_Res
      CPI.DT_Res$CPI.DT_Avg <- apply(CPI.DT_Res_orig, 1, mean, na.rm=TRUE)
      CPI.DT_Res$CPI.DT_AvgEx8 <- apply(CPI.DT_Res_orig %>% dplyr::select(-ends_with("j8")), 1, mean, na.rm=TRUE)
      for (i in 1:6) {
        CPI.DT_Res[,paste0("CPI.DT_D",i)] <- apply(CPI.DT_Res_orig %>% dplyr::select(contains(paste0("a",i))), 
                                                   1, mean, na.rm=TRUE)
      }
      for (i in 2:7) {
        CPI.DT_Res[,paste0("CPI.DT_T",i)] <- apply(CPI.DT_Res_orig %>% dplyr::select(contains(paste0("j",i))), 
                                                   1, mean, na.rm=TRUE)
      }
      CPI.DT_Res$CPI.DT_DAvg <- apply(CPI.DT_Res %>% dplyr::select(starts_with("CPI.DT_D")),
                                      1, mean, na.rm=TRUE)
      CPI.DT_Res$CPI.DT_TAvg <- apply(CPI.DT_Res %>% dplyr::select(starts_with("CPI.DT_T")),
                                      1, mean, na.rm=TRUE)
      Res.perm <- cbind(Res.perm, as.matrix(CPI.DT_Res))
    }
    if (GEE) {
      GEE_Res <- sapply(Perm.Dat.Long,
                        FUN=function(x) unname(geeglm(Y~Interv+factor(Period),
                                                      data=x %>% dplyr::mutate(Clf=factor(Cluster)),
                                                      id=Clf,
                                                      corstr=corstr)$coefficients["Interv"]))
      Res.perm <- cbind(Res.perm, Comp_GEE=GEE_Res)
    }
  }
  return(Res.perm)
}

### Framework for simulated data with fixed treatment schedule:
simulate_SWT <- function(NumSims,
                          N, J, StartingPds=NULL,
                          mu, Alpha1, 
                          T1, T2, ProbT1,
                          sig_nu, sig_e, m,
                         ThetaType, ThetaDF,
                         MVO_list, SO_list=NULL,
                         Comparisons=NULL, corstr="exchangeable",
                         Permutations=0) {
  if (NumSims <= 1) {
    stop(simpleError(message="NumSims must be at least 2"))
  }
  Sim.Fr <- Sim_Frame(N, J, StartingPds)
  Sim.Dat <- replicate(n=NumSims, 
                       as.matrix(Sim_Data(Sim.Fr, mu, Alpha1, 
                                          T1, T2, ProbT1,
                                          sig_nu, sig_e, m,
                                          ThetaType, ThetaDF)), 
                       simplify="array")
  Sim.Wt <- Sim_Weights(MV_Out=MVO_list,
                        Solve_Out=SO_list, 
                        Comparisons=Comparisons)
  Sim.Res <- Sim_Analyze(Sim.Dat, Sim.Wt,
                         MEM=("MEM" %in% Comparisons),
                         CPI=("CPI" %in% Comparisons),
                         CPI.T=("CPI.T" %in% Comparisons),
                         CPI.D=("CPI.D" %in% Comparisons),
                         CPI.DT=("CPI.DT" %in% Comparisons),
                         GEE=("GEE" %in% Comparisons),
                         corstr=corstr)
  if (Permutations > 0) {
    Sim.Perm.Res <- replicate(n=Permutations,
                              expr=Sim_Permutation(Sim.Dat, Sim.Wt,
                                    N, J,
                                    MEM=("MEM" %in% Comparisons),
                                    CPI=("CPI" %in% Comparisons),
                                    CPI.T=("CPI.T" %in% Comparisons),
                                    CPI.D=("CPI.D" %in% Comparisons),
                                    CPI.DT=("CPI.DT" %in% Comparisons),
                                    GEE=("GEE" %in% Comparisons),
                                    corstr=corstr),
                              simplify="array")
    Sim.All <- array(c(Sim.Res, Sim.Perm.Res),
                     dim=c(dim(Sim.Res)[1:2],dim(Sim.Perm.Res)[3]+1))
    dimnames(Sim.All) <- list(paste("Simulation",1:(dim(Sim.All)[1]),sep="_"),
                              colnames(Sim.Res),
                              c("Results",paste("Placebo",1:(dim(Sim.Perm.Res)[3]),sep="_")))
    Sim.PVals <- apply(X=Sim.All, MARGIN=c(1,2),
                       FUN=function(x) mean(abs(x["Results"]) <= abs(x), na.rm=TRUE))
    return(list(Estimates=Sim.All[,,"Results"],
                PValues=Sim.PVals))
  } else {
    dimnames(Sim.Res) <- list(paste("Simulation",1:(dim(Sim.Res)[1]),sep="_"),
                              colnames(Sim.Res))
    return(list(Estimates=Sim.Res))
  }
}

simulate_FromSet <- function(Param_Set,
                             Theta_Set,
                             StartingPds=NULL,
                             Alpha1, 
                             T1, T2, 
                             MVO_list, SO_list=NULL,
                             outdir=NULL,
                             outname=NULL) {
      row <- Param_Set[i,]
      print(paste("Starting Sim. Number",row$SimNo))
      assign(x=paste0("Res_Sim_",i),
             value=simulate_SWT(row$NumSims,
                                row$N, row$J, StartingPds,
                                row$mu, Alpha1, 
                                T1, T2, row$ProbT1,
                                row$sig_nu, row$sig_e, row$m,
                                Theta_Set[[i]]$Type, Theta_Set[[i]]$ThetaDF,
                                MVO_list, SO_list,
                                Comparisons=Theta_Set[[i]]$Comps, Theta_Set[[i]]$corstr,
                                Permutations=row$NumPerms))
      save(list=paste0("Res_Sim_",row$SimNo), 
           file=paste0(outdir,"/",outname,"_",row$SimNo,".Rda"))
}
