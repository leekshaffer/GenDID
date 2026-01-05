#######################################
###### File: Sim_Single.R #############
###### Lee Kennedy-Shaffer ############
#######################################

source("R/Full_Analysis.R")
source("R/CompEsts.R")

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

  ### Permuted Orders:
  P.Orders <- NULL
  if (NumPerms > 0) {
    P.Orders <- replicate(n=NumPerms,
                          expr=Permute_order(N, J),
                          simplify=FALSE)
  }

  ## Get GenDID and Nested Method Estimates and Inference

  ### Full list of observation weights for assumption settings
  Obs.W <- do.call("cbind",
          lapply(1:length(MVO_list),
                 function(x) as_tibble(MVO_list[[x]]$MV$Obs.weights) %>%
                   dplyr::rename_with(~paste(names(MVO_list)[x],
                                             colnames(MVO_list[[x]]$MV$Obs.weights),
                                             sep="_",recycle0=FALSE))))
  if (!is.null(Comps.Nest)) {
    Obs.W <- cbind(Obs.W, Comp_wts)
  }

  ### Estimates calculation:
  Estimates <- t(Obs.W) %*% Obs_Y

  ### P Values calculation:
  if (!is.null(P.Orders)) {
    Perm_Ests <- lapply(P.Orders, FUN=function(x) t(Obs.W) %*% Obs_Y[x,])
    PVals <- apply(simplify2array(
      lapply(Perm_Ests,
             FUN=function(x) abs(x) >= abs(Estimates))),
      MARGIN=c(1,2), mean)
  }

  ### CIs calculation:
  if (!is.null(P.Orders) & !is.null(CI.GDID)) {
    CI_Ps <- lapply(names(CI.GDID),
                    FUN=function(CI.name) {
                      CI.FX <- as.matrix(CI.GDID[[CI.name]])
                      CI.Estimates <- Estimates - t(Obs.W) %*% CI.FX
                      CI_Perm_Ests <- lapply(P.Orders,
                                             FUN=function(x) t(Obs.W) %*% (Obs_Y - CI.FX)[x,])
                      CI_Ps <- apply(simplify2array(
                        lapply(CI_Perm_Ests,
                               FUN=function(x) abs(x) >= abs(CI.Estimates))),
                        MARGIN=c(1,2), mean)
                    })
    names(CI_Ps) <- names(CI.GDID)
    if (is.null(CI.SV.Mesh)) {
      CI_List <- lapply(CI_Ps, FUN=function(x) x >= 1-CI.Perc)
      CI_Mesh <- NULL
    } else {
      CI_List <- lapply(CI_Ps[names(CI_Ps)[!(names(CI_Ps) %in% CI.SV.Mesh$Name)]],
                        FUN=function(x) x >= 1-CI.Perc)
      CI_Results <- do.call("rbind",
                            lapply(CI.SV.Mesh$Name,
                                   FUN=function(x) as_tibble(CI_Ps[[x]], rownames="Estimator") %>%
                                     dplyr::mutate(Name=x, EstNum=row_number()))) %>%
        pivot_longer(cols=-c(EstNum,Estimator,Name), names_to="Outcome", values_to="CI_P") %>%
        left_join(CI.SV.Mesh %>% pivot_longer(cols=-c(ValNum, Name), names_to="Outcome", values_to="Effect"),
                  by=join_by(Name,Outcome)) %>%
        dplyr::group_by(EstNum,Estimator,Outcome) %>%
        dplyr::filter(CI_P >= 1-CI.Perc) %>%
        dplyr::summarize(CIL=min(Effect),
                         CIU=max(Effect),
                         .groups="drop") %>%
        dplyr::mutate(CIW=CIU-CIL)
      CI_Mesh <- list(CIL=CI_Results %>% dplyr::select(Estimator,Outcome,CIL) %>%
                        pivot_wider(id_cols=Estimator, names_from=Outcome, values_from=CIL) %>%
                        dplyr::select(all_of(c("Estimator",colnames(PVals)))),
                      CIU=CI_Results %>% dplyr::select(Estimator,Outcome,CIU) %>%
                        pivot_wider(id_cols=Estimator, names_from=Outcome, values_from=CIU) %>%
                        dplyr::select(all_of(c("Estimator",colnames(PVals)))),
                      CIW=CI_Results %>% dplyr::select(Estimator,Outcome,CIW) %>%
                        pivot_wider(id_cols=Estimator, names_from=Outcome, values_from=CIW) %>%
                        dplyr::select(all_of(c("Estimator",colnames(PVals)))))

    }
  }

  ### Comparisons

  Data.ForLME <- Data %>%
    dplyr::select(Cluster,Period,Start,Interv,all_of(starts_with("Y.ij."))) %>%
    dplyr::select(-c(Y.ij.bar,Y.ij.sd)) %>%
    dplyr::mutate(PeriodF=factor(Period),
                  ClusterF=factor(Cluster),
                  CPI=factor(paste(Cluster,Period,sep="_")),
                  Diff=if_else(Period-Start < 0, 0, Period-Start+1),
                  DiffF=factor(Diff))

  Data.Long <- Data.ForLME %>%
    tidyr::pivot_longer(cols=starts_with("Y.ij."),
                        names_to="Indiv", names_prefix="Y.ij.",
                        values_to="Outcome")


  Data.Yijs <- Data.ForLME %>% dplyr::select(all_of(starts_with("Y.ij.")))

  Data.Fr <- Data.ForLME %>%
    dplyr::select(Cluster,Period,Start,Interv,PeriodF,ClusterF,CPI,Diff,DiffF)


  if (!is.null(Comps_PermPs)) {
    Data.Perms <- lapply(P.Orders,
                         FUN=function(x) Data.Fr %>% bind_cols(Data.Yijs[x,]) %>%
                           tidyr::pivot_longer(cols=starts_with("Y.ij."),
                                               names_to="Indiv", names_prefix="Y.ij.",
                                               values_to="Outcome"))
    Data.Perms.Short <- lapply(P.Orders,
                               FUN=function(x) Data.Fr %>% bind_cols(Data[x,"Y.ij.bar"]))
  }

  Comp_Outs <- NULL

  if ("TW" %in% Comps) {
    TW <- summary(feols(Y.ij.bar~Interv | Cluster + Period,
                        data=Data),
                  vcov=~Cluster)
    TW.row <- tibble_row(Method="TW",
                         Estimate=TW$coeftable["Interv","Estimate"],
                         SE=TW$coeftable["Interv","Std. Error"],
                         P=TW$coeftable["Interv","Pr(>|t|)"],
                         CIL=confint(TW, level=CI.Perc)["Interv",1],
                         CIU=confint(TW, level=CI.Perc)["Interv",2])
    if ("TW" %in% Comps_PermPs) {
      if ("W_TW" %in% rownames(PVals)) {
        TW.row <- TW.row %>% dplyr::mutate(P.Perm=PVals["W_TW",1])
      } else {
        TW.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) summary(feols(Y.ij.bar~Interv | Cluster + Period,
                                                                        data=x),
                                                                  vcov=~Cluster)$coeftable["Interv","Estimate"]))
        TW.row <- TW.row %>% dplyr::mutate(P.Perm=mean(abs(TW.Perms) >= abs(TW.row$Estimate)))
      }
    }
    Comp_Outs <- Comp_Outs %>%
      bind_rows(TW.row)
  }

  if ("CS" %in% Comps) {
    CS_gt <- att_gt(yname="Y.ij.bar",
                    tname="Period",
                    idname="Cluster",
                    gname="Start",
                    data=Data,
                    panel=TRUE,
                    control_group="notyettreated")
    simple <- aggte(CS_gt, type="simple")
    dynamic <- aggte(CS_gt, type="dynamic")
    group <- aggte(CS_gt, type="group")
    cal <- aggte(CS_gt, type="calendar")

    CS.rows <- tibble(Method=paste0("CS_",c("Simple","Dynamic","Group","Calendar")),
                      Estimate=c(simple$overall.att,
                                 dynamic$overall.att,
                                 group$overall.att,
                                 cal$overall.att),
                      SE=c(simple$overall.se,
                           dynamic$overall.se,
                           group$overall.se,
                           cal$overall.se)) %>%
      dplyr::mutate(P=2*pnorm(abs(Estimate/SE), lower.tail=FALSE),
                    CIL=Estimate+qnorm((1-CI.Perc)/2)*SE,
                    CIU=Estimate+qnorm((1+CI.Perc)/2)*SE)

    if ("CS" %in% Comps_PermPs) {
      P.Perms <- NULL
      if ("W_CS.W_simple" %in% rownames(PVals)) {
        P.Perms <- c(P.Perms, PVals["W_CS.W_simple",1])
      } else {
        CS.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) aggte(att_gt(yname="Y.ij.bar",
                                                                          tname="Period",
                                                                          idname="Cluster",
                                                                          gname="Start",
                                                                          data=x,
                                                                          panel=TRUE,
                                                                          control_group="notyettreated"),
                                                                type="simple")$overall.att))
        P.Perms <- c(P.Perms, mean(abs(CS.Perms) >= abs(CS.rows %>% dplyr::filter(Method=="CS_Simple") %>% pull(Estimate))))
      }

      if ("W_CS.W_dynamic" %in% rownames(PVals)) {
        P.Perms <- c(P.Perms, PVals["W_CS.W_dynamic",1])
      } else {
        CS.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) aggte(att_gt(yname="Y.ij.bar",
                                                                       tname="Period",
                                                                       idname="Cluster",
                                                                       gname="Start",
                                                                       data=x,
                                                                       panel=TRUE,
                                                                       control_group="notyettreated"),
                                                                type="dynamic")$overall.att))
        P.Perms <- c(P.Perms, mean(abs(CS.Perms) >= abs(CS.rows %>% dplyr::filter(Method=="CS_Dynamic") %>% pull(Estimate))))
      }
    }

    if ("W_CS.W_group" %in% rownames(PVals)) {
      P.Perms <- c(P.Perms, PVals["W_CS.W_group",1])
    } else {
      CS.Perms <- simplify2array(lapply(Data.Perms.Short,
                                        FUN=function(x) aggte(att_gt(yname="Y.ij.bar",
                                                                     tname="Period",
                                                                     idname="Cluster",
                                                                     gname="Start",
                                                                     data=x,
                                                                     panel=TRUE,
                                                                     control_group="notyettreated"),
                                                              type="group")$overall.att))
      P.Perms <- c(P.Perms, mean(abs(CS.Perms) >= abs(CS.rows %>% dplyr::filter(Method=="CS_Group") %>% pull(Estimate))))
    }

    if ("W_CS.W_calendar" %in% rownames(PVals)) {
      P.Perms <- c(P.Perms, PVals["W_CS.W_calendar",1])
    } else {
      CS.Perms <- simplify2array(lapply(Data.Perms.Short,
                                        FUN=function(x) aggte(att_gt(yname="Y.ij.bar",
                                                                     tname="Period",
                                                                     idname="Cluster",
                                                                     gname="Start",
                                                                     data=x,
                                                                     panel=TRUE,
                                                                     control_group="notyettreated"),
                                                              type="calendar")$overall.att))
      P.Perms <- c(P.Perms, mean(abs(CS.Perms) >= abs(CS.rows %>% dplyr::filter(Method=="CS_Calendar") %>% pull(Estimate))))
    }

    Comp_Outs <- Comp_Outs %>%
      bind_rows(CS.rows %>% bind_cols(P.Perm=P.Perms))
  }

  if ("SA" %in% Comps) {
    ## Data without the max period to create a never-treated group:
    Max_Pd <- max(Data$Period)
    Data.SA <- Data %>% filter(Period != Max_Pd) %>%
      mutate(Start=if_else(Start==Max_Pd, 0, Start))
    SA <- feols(Y.ij.bar~sunab(cohort=Start,
                              period=Period,
                              ref.c=NULL,
                              ref.p=c(-1,Max_Pd-2), att=TRUE) | Cluster + Period,
                data=Data.SA)
    SA.CI <- confint(SA, level=CI.Perc)
    SA.row <- tibble_row(Method="SA",
                         Estimate=SA$coeftable["ATT","Estimate"],
                         SE=SA$coeftable["ATT","Std. Error"],
                         P=SA$coeftable["ATT","Pr(>|t|)"],
                         CIL=SA.CI["ATT",1],
                         CIU=SA.CI["ATT",2])

    if ("SA" %in% Comps_PermPs) {
      if ("W_SA.W_ATT" %in% rownames(PVals)) {
        SA.row <- SA.row %>% dplyr::mutate(P.Perm=PVals["W_SA.W_ATT",1])
      } else {
        SA.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) feols(Y.ij.bar~sunab(cohort=Start,
                                                                                period=Period,
                                                                                ref.c=NULL,
                                                                                ref.p=c(-1,Max_Pd-2), att=TRUE) | Cluster + Period,
                                                                 data = x %>% filter(Period != Max_Pd) %>%
                                                                  mutate(Start=if_else(Start==Max_Pd, 0, Start)))$coeftable["ATT","Estimate"]))
        SA.row <- SA.row %>% dplyr::mutate(P.Perm=mean(abs(SA.Perms) >= abs(SA.row$Estimate)))
      }
    }
    Comp_Outs <- Comp_Outs %>%
      bind_rows(SA.row)
  }

  if ("CH" %in% Comps) {
    CH <- did_multiplegt(
      mode="dyn",
      df=Data,
      outcome="Y.ij.bar",
      group="Cluster",
      time="Period",
      treatment="Interv",
      graph_off=TRUE,
      ci_level=CI.Perc*100
    )$results$Effects

    CH.row <- tibble_row(Method="CH",
                         Estimate=CH[1,"Estimate"],
                         SE=CH[1,"SE"],
                         CIL=CH[1,"LB CI"],
                         CIU=CH[1,"UB CI"]) %>%
      dplyr::mutate(P=2*pnorm(abs(Estimate/SE), lower.tail=FALSE))

    if ("CH" %in% Comps_PermPs) {
      if ("W_CH.W_M" %in% rownames(PVals)) {
        CH.row <- CH.row %>% dplyr::mutate(P.Perm=PVals["W_CH.W_M",1])
      } else {
        CH.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) (did_multiplegt(
                                            mode="dyn",
                                            df=Data,
                                            outcome="Y.ij.bar",
                                            group="Cluster",
                                            time="Period",
                                            treatment="Interv",
                                            graph_off=TRUE,
                                            ci_level=CI.Perc*100
                                          )$results$Effects)[1,"Estimate"]))
        CH.row <- CH.row %>% dplyr::mutate(P.Perm=mean(abs(CH.Perms) >= abs(CH.row$Estimate)))
      }
    }
    Comp_Outs <- Comp_Outs %>%
      bind_rows(CH.row)
  }

  if ("MEM" %in% Comps) {
    MEM <- lmer(Outcome~Interv+PeriodF+(1|ClusterF), data=Data.Long)
    MEM.ci <- confint(MEM, parm="Interv", level=CI.Perc, quiet=TRUE)
    MEM.row <- tibble_row(Method="MEM",
                          Estimate=coeftable(MEM)["Interv","Estimate"],
                          SE=coeftable(MEM)["Interv","Std. Error"],
                          CIL=MEM.ci["Interv",1],
                          CIU=MEM.ci["Interv",2]) %>%
      dplyr::mutate(P=2*pnorm(abs(Estimate/SE), lower.tail=FALSE))
    if ("MEM" %in% Comps_PermPs) {
      MEM.Perms <- simplify2array(lapply(Data.Perms,
                          FUN=function(x) coeftable(lmer(Outcome~Interv+PeriodF+(1|ClusterF), data=x))["Interv","Estimate"]))
      MEM.row <- MEM.row %>% dplyr::mutate(P.Perm=mean(abs(MEM.Perms) >= abs(MEM.row$Estimate)))
    }
    Comp_Outs <- Comp_Outs %>% bind_rows(MEM.row)
  }

  if ("CPI" %in% Comps) {
    CPI <- lmer(Outcome~Interv+PeriodF+(1|ClusterF)+(1|CPI), data=Data.Long,
                control=lmerControl(check.rankX="silent.drop.cols",
                                    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    CPI.ci <- confint(CPI, parm="Interv", level=CI.Perc, quiet=TRUE)
    CPI.row <- tibble_row(Method="CPI",
                          Estimate=coeftable(CPI)["Interv","Estimate"],
                          SE=coeftable(CPI)["Interv","Std. Error"],
                          CIL=CPI.ci["Interv",1],
                          CIU=CPI.ci["Interv",2]) %>%
      dplyr::mutate(P=2*pnorm(abs(Estimate/SE), lower.tail=FALSE))
    if ("CPI" %in% Comps_PermPs) {
      CPI.Perms <- simplify2array(lapply(Data.Perms,
                                         FUN=function(x) coeftable(lmer(Outcome~Interv+PeriodF+(1|ClusterF)+(1|CPI), data=x,
                                                                        control=lmerControl(check.rankX="silent.drop.cols",
                                                                                            check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))["Interv","Estimate"]))
      CPI.row <- CPI.row %>% dplyr::mutate(P.Perm=mean(abs(CPI.Perms) >= abs(CPI.row$Estimate)))
    }
    Comp_Outs <- Comp_Outs %>% bind_rows(CPI.row)
  }

  if ("CPI.T" %in% Comps) {
    CPI.T <- lmer(Outcome~Interv+Interv:PeriodF+PeriodF+(1|ClusterF)+(1|CPI), data=Data.Long,
                  control=lmerControl(check.rankX="silent.drop.cols",
                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    IntVec <- grepl("Interv:",rownames(coeftable(CPI.T)))
    IVec <- (rownames(coeftable(CPI.T))=="Interv")+IntVec/(sum(IntVec)+1)
    CPI.T.row <- tibble_row(Method="CPI.T",
                            Estimate=(IVec %*% coeftable(CPI.T)[,"Estimate"])[1,1],
                            SE=sqrt(IVec %*% vcov(CPI.T) %*% IVec)[1,1]) %>%
      dplyr::mutate(P=2*pnorm(abs(Estimate/SE), lower.tail=FALSE),
                    CIL=Estimate+qnorm((1-CI.Perc)/2)*SE,
                    CIU=Estimate+qnorm((1+CI.Perc)/2)*SE)
    if ("CPI.T" %in% Comps_PermPs) {
      CPI.T.Perms <- simplify2array(lapply(Data.Perms,
                                           FUN=function(x) (IVec %*% coeftable(lmer(Outcome~Interv+Interv:PeriodF+PeriodF+(1|ClusterF)+(1|CPI), data=x,
                                                                                    control=lmerControl(check.rankX="silent.drop.cols",
                                                                                                        check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))[,"Estimate"])[1,1]))
      CPI.T.row <- CPI.T.row %>% dplyr::mutate(P.Perm=mean(abs(CPI.T.Perms) >= abs(CPI.T.row$Estimate)))
    }
    Comp_Outs <- Comp_Outs %>% bind_rows(CPI.T.row)
  }

  if ("CPI.D" %in% Comps) {
    CPI.D <- lmer(Outcome~Interv+Interv:DiffF+PeriodF+(1|ClusterF)+(1|CPI), data=Data.Long,
                control=lmerControl(check.rankX="silent.drop.cols",
                                    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    IntVec <- grepl("Interv:",rownames(coeftable(CPI.D)))
    IVec <- (rownames(coeftable(CPI.D))=="Interv")+IntVec/(sum(IntVec)+1)
    CPI.D.row <- tibble_row(Method="CPI.D",
                          Estimate=(IVec %*% coeftable(CPI.D)[,"Estimate"])[1,1],
                          SE=sqrt(IVec %*% vcov(CPI.D) %*% IVec)[1,1]) %>%
      dplyr::mutate(P=2*pnorm(abs(Estimate/SE), lower.tail=FALSE),
                    CIL=Estimate+qnorm((1-CI.Perc)/2)*SE,
                    CIU=Estimate+qnorm((1+CI.Perc)/2)*SE)
    if ("CPI.D" %in% Comps_PermPs) {
      CPI.D.Perms <- simplify2array(lapply(Data.Perms,
                                         FUN=function(x) (IVec %*% coeftable(lmer(Outcome~Interv+Interv:DiffF+PeriodF+(1|ClusterF)+(1|CPI), data=x,
                                                                        control=lmerControl(check.rankX="silent.drop.cols",
                                                                                            check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))[,"Estimate"])[1,1]))
      CPI.D.row <- CPI.D.row %>% dplyr::mutate(P.Perm=mean(abs(CPI.D.Perms) >= abs(CPI.D.row$Estimate)))
    }
    Comp_Outs <- Comp_Outs %>% bind_rows(CPI.D.row)
  }

  if ("CPI.DT" %in% Comps) {
    CPI.DT <- lmer(Outcome~Interv+Interv:DiffF:PeriodF+PeriodF+(1|ClusterF)+(1|CPI), data=Data.Long,
                  control=lmerControl(check.rankX="silent.drop.cols",
                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    IntVec <- grepl("Interv:",rownames(coeftable(CPI.DT)))
    IVec <- (rownames(coeftable(CPI.DT))=="Interv")+IntVec/(sum(IntVec)+1)
    CPI.DT.row <- tibble_row(Method="CPI.DT",
                            Estimate=(IVec %*% coeftable(CPI.DT)[,"Estimate"])[1,1],
                            SE=sqrt(IVec %*% vcov(CPI.DT) %*% IVec)[1,1]) %>%
      dplyr::mutate(P=2*pnorm(abs(Estimate/SE), lower.tail=FALSE),
                    CIL=Estimate+qnorm((1-CI.Perc)/2)*SE,
                    CIU=Estimate+qnorm((1+CI.Perc)/2)*SE)
    if ("CPI.DT" %in% Comps_PermPs) {
      CPI.DT.Perms <- simplify2array(lapply(Data.Perms,
                                           FUN=function(x) (IVec %*% coeftable(lmer(Outcome~Interv+Interv:DiffF:PeriodF+PeriodF+(1|ClusterF)+(1|CPI), data=x,
                                                                                    control=lmerControl(check.rankX="silent.drop.cols",
                                                                                                        check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))[,"Estimate"])[1,1]))
      CPI.DT.row <- CPI.DT.row %>% dplyr::mutate(P.Perm=mean(abs(CPI.DT.Perms) >= abs(CPI.DT.row$Estimate)))
    }
    Comp_Outs <- Comp_Outs %>% bind_rows(CPI.DT.row)
  }

  if (("CLWP" %in% Comps) | ("CLWPA" %in% Comps)) {
    if ("CPI" %in% Comps) {
      start_theta <- CPI.row$Estimate
      start_sigma <- sqrt(as.data.frame(VarCorr(CPI))[2,"vcov"]*2)
    } else if ("MEM" %in% Comps) {
      start_theta <- MEM.row$Estimate
      start_sigma <- sqrt(as.data.frame(VarCorr(MEM))[2,"vcov"]*2)
    } else {
      MEM <- lmer(Outcome~Interv+PeriodF+(1|ClusterF), data=Data.Long)
      start_theta <- coeftable(MEM)["Interv","Estimate"]
      start_sigma <- sqrt(as.data.frame(VarCorr(MEM))[2,"vcov"]*2)
    }
    if ("CLWP" %in% Comps) {
      CLWP <- CLWP_fit(my.data=Data,
                       start_theta=start_theta,
                       start_sigma=start_sigma,
                       N=N)
      CLWP.row <- tibble_row(Method="CLWP",
                             Estimate=CLWP["CLWP_est.theta"],
                             SE=CLWP["CLWP_se"],
                             P=CLWP["CLWP_pval.theta"]) %>%
        dplyr::mutate(CIL=Estimate+qnorm((1-CI.Perc)/2)*SE,
                      CIU=Estimate+qnorm((1+CI.Perc)/2)*SE)
      if ("CLWP" %in% Comps_PermPs) {
        CLWP.Perms <- simplify2array(lapply(Data.Perms.Short,
                                              FUN=function(x) unname(CLWP_fit(my.data=x,
                                                                       start_theta=start_theta,
                                                                       start_sigma=start_sigma,
                                                                       N=N)["CLWP_est.theta"])))
        CLWP.row <- CLWP.row %>% dplyr::mutate(P.Perm=mean(abs(CLWP.Perms) >= abs(CLWP.row$Estimate)))
      }
      Comp_Outs <- Comp_Outs %>% bind_rows(CLWP.row)
    }

    if ("CLWPA" %in% Comps) {
      CLWPA <- CLWPA_fit(my.data=Data,
                       start_theta=start_theta,
                       start_sigma=start_sigma,
                       N=N)
      CLWPA.row <- tibble_row(Method="CLWPA",
                             Estimate=CLWPA["CLWPA_est.theta"],
                             SE=CLWPA["CLWPA_se"],
                             P=CLWPA["CLWPA_pval.theta"]) %>%
        dplyr::mutate(CIL=Estimate+qnorm((1-CI.Perc)/2)*SE,
                      CIU=Estimate+qnorm((1+CI.Perc)/2)*SE)
      if ("CLWPA" %in% Comps_PermPs) {
        CLWPA.Perms <- simplify2array(lapply(Data.Perms.Short,
                                            FUN=function(x) unname(CLWPA_fit(my.data=x,
                                                                            start_theta=start_theta,
                                                                            start_sigma=start_sigma,
                                                                            N=N)["CLWPA_est.theta"])))
        CLWPA.row <- CLWPA.row %>% dplyr::mutate(P.Perm=mean(abs(CLWPA.Perms) >= abs(CLWPA.row$Estimate)))
      }
      Comp_Outs <- Comp_Outs %>% bind_rows(CLWPA.row)
    }
  }

  return(c(list(Estimates=Estimates,
                P_Values=PVals),
           CI_List=CI_List,
           CI_Mesh=CI_Mesh,
           list(Comparisons=Comp_Outs)))
}

# Speed Testing:

## From test on Work Comp: 23 seconds for 20 permutations w/ mesh of 25
## From test on Work Comp: 28 seconds for 20 permutations w/ mesh of 501
## From test on Work Comp: 175 seconds for 100 permutations w/ mesh of 501
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
