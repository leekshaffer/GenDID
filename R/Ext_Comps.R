#######################################
###### File: Ext_Comps.R ##############
#######################################

require(dplyr)
require(tidyr)
require(tibble)
require(lme4) ## For mixed effects models and CLWP
require(bbmle) ## For CLWP

source("R/CLWP_Fns.R")
source("R/Analysis.R") ## For permutation function

## Need long data with Outcome as the individual-level outcome column
## Note: SummOutName specifies an existing cluster-period summary column to use;
## Otherwise the average of the Outcome column is used as a summary
## Also need Cluster, Period, Start, Interv columns

Ext_Comps <- function(Data.Long,
                      SummOutName=NULL,
                      Comps=c("TW","CS","SA","CH","MEM",
                              "CPI","CPI.T","CPI.D","CPI.DT",
                              "CLWP","CLWPA"),
                      Comps_PermPs=c("TW","CS","SA","CH","MEM",
                                     "CPI","CPI.T","CPI.D","CPI.DT",
                                     "CLWP","CLWPA"),
                      P.Orders=NULL,
                      NumPerms=0,
                      Results=NULL,
                      CI.Perc=0.95) {

  Data.LME <- Data.Long %>%
    dplyr::mutate(PeriodF=factor(Period),
                  ClusterF=factor(Cluster),
                  CPI=factor(paste(Cluster,Period,sep="_")),
                  Diff=if_else(Period-Start < 0, 0, Period-Start+1),
                  DiffF=factor(Diff))

  if (!is.null(SummOutName)) {
    Data <- Data.LME %>%
      dplyr::group_by(Cluster,Period,Start,Interv,
                      PeriodF,ClusterF,CPI,Diff,DiffF) %>%
      dplyr::summarize(Summ=mean(.data[[SummOutName]], na.rm=TRUE),
                       .groups="drop") %>%
      dplyr::mutate(CPN=row_number())
  } else {
    Data <- Data.LME %>%
      dplyr::group_by(Cluster,Period,Start,Interv,
                      PeriodF,ClusterF,CPI,Diff,DiffF) %>%
      dplyr::summarize(Summ=mean(Outcome, na.rm=TRUE),
                       .groups="drop") %>%
      dplyr::mutate(CPN=row_number())
  }

  N <- length(unique(Data$Cluster))
  J <- length(unique(Data$Period))

  Data.LME <- Data.LME %>%
    left_join(Data %>%
                dplyr::select(Cluster,Period,CPN),
              by=join_by(Cluster,Period))

  Data.Fr <- Data %>%
    dplyr::select(-Summ)

  Data.Summ <- Data %>%
    dplyr::select(CPN,Summ)

  Data.Outcome <- Data.LME %>%
    dplyr::select(CPN,Outcome)

  if (!is.null(Comps_PermPs)) {
    if (is.null(P.Orders)) {
      if (NumPerms > 0) {
        P.Orders <- replicate(n=NumPerms,
                              expr=Permute_order(N, J),
                              simplify=FALSE)
      } else {
        P.Orders <- NULL
      }
    }
    Data.Perms <- lapply(P.Orders,
                         FUN=function(x) {
                           Data.Outcome %>%
                             dplyr::left_join(tibble(CPN=x) %>% dplyr::mutate(NewCPN=row_number()),
                                              by=join_by(CPN)) %>%
                             dplyr::select(-CPN) %>%
                             dplyr::rename(CPN=NewCPN) %>%
                             dplyr::left_join(Data.Fr, by=join_by(CPN)) %>%
                             dplyr::arrange(CPN)
                         })
    Data.Perms.Short <- lapply(P.Orders,
                               FUN=function(x) {
                                 Data.Summ %>%
                                   dplyr::left_join(tibble(CPN=x) %>% dplyr::mutate(NewCPN=row_number()),
                                                    by=join_by(CPN)) %>%
                                   dplyr::select(-CPN) %>%
                                   dplyr::rename(CPN=NewCPN) %>%
                                   dplyr::left_join(Data.Fr, by=join_by(CPN)) %>%
                                   dplyr::arrange(CPN)
                               })
  }

  Comp_Outs <- NULL

  if ("TW" %in% Comps) {
    require(fixest)
    TW <- summary(feols(Summ~Interv | Cluster + Period,
                        data=Data),
                  vcov=~Cluster)
    TW.row <- tibble_row(Method="TW",
                         Estimate=TW$coeftable["Interv","Estimate"],
                         SE=TW$coeftable["Interv","Std. Error"],
                         P=TW$coeftable["Interv","Pr(>|t|)"],
                         CIL=confint(TW, level=CI.Perc)["Interv",1],
                         CIU=confint(TW, level=CI.Perc)["Interv",2])
    if ("TW" %in% Comps_PermPs) {
      if ("W_TW" %in% Results$Estimator) {
        TW.row <- TW.row %>% dplyr::mutate(P.Perm=Results %>%
                                             dplyr::filter(Estimator=="W_TW",Outcome=="Probability") %>%
                                             pull(P))
      } else {
        TW.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) summary(feols(Summ~Interv | Cluster + Period,
                                                                        data=x),
                                                                  vcov=~Cluster)$coeftable["Interv","Estimate"]))
        TW.row <- TW.row %>% dplyr::mutate(P.Perm=mean(abs(TW.Perms) >= abs(TW.row$Estimate)))
      }
    }
    Comp_Outs <- Comp_Outs %>%
      bind_rows(TW.row)
  }

  if ("CS" %in% Comps) {
    require(did)
    CS_gt <- att_gt(yname="Summ",
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
      if (sum(c("W_CS.W_simple","W_CS.W_dynamic","W_CS.W_group","W_CS.W_calendar") %in% Results$Estimator)==4) {
        P.Perms <- c(Results %>%
                       dplyr::filter(Estimator=="W_CS.W_simple",
                                     Outcome=="Probability") %>%
                       pull(P),
                     Results %>%
                       dplyr::filter(Estimator=="W_CS.W_dynamic",
                                     Outcome=="Probability") %>%
                       pull(P),
                     Results %>%
                       dplyr::filter(Estimator=="W_CS.W_group",
                                     Outcome=="Probability") %>%
                       pull(P),
                     Results %>%
                       dplyr::filter(Estimator=="W_CS.W_calendar",
                                     Outcome=="Probability") %>%
                       pull(P))
      } else {
        P.Perms <- NULL
        CS.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) aggte(att_gt(yname="Summ",
                                                                       tname="Period",
                                                                       idname="Cluster",
                                                                       gname="Start",
                                                                       data=x,
                                                                       panel=TRUE,
                                                                       control_group="notyettreated"),
                                                                type="simple")$overall.att))
        P.Perms <- c(P.Perms, mean(abs(CS.Perms) >= abs(CS.rows %>% dplyr::filter(Method=="CS_Simple") %>% pull(Estimate))))
        CS.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) aggte(att_gt(yname="Summ",
                                                                       tname="Period",
                                                                       idname="Cluster",
                                                                       gname="Start",
                                                                       data=x,
                                                                       panel=TRUE,
                                                                       control_group="notyettreated"),
                                                                type="dynamic")$overall.att))
        P.Perms <- c(P.Perms, mean(abs(CS.Perms) >= abs(CS.rows %>% dplyr::filter(Method=="CS_Dynamic") %>% pull(Estimate))))
        CS.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) aggte(att_gt(yname="Summ",
                                                                       tname="Period",
                                                                       idname="Cluster",
                                                                       gname="Start",
                                                                       data=x,
                                                                       panel=TRUE,
                                                                       control_group="notyettreated"),
                                                                type="group")$overall.att))
        P.Perms <- c(P.Perms, mean(abs(CS.Perms) >= abs(CS.rows %>% dplyr::filter(Method=="CS_Group") %>% pull(Estimate))))
        CS.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) aggte(att_gt(yname="Summ",
                                                                       tname="Period",
                                                                       idname="Cluster",
                                                                       gname="Start",
                                                                       data=x,
                                                                       panel=TRUE,
                                                                       control_group="notyettreated"),
                                                                type="calendar")$overall.att))
        P.Perms <- c(P.Perms, mean(abs(CS.Perms) >= abs(CS.rows %>% dplyr::filter(Method=="CS_Calendar") %>% pull(Estimate))))
      }
    }

    Comp_Outs <- Comp_Outs %>%
      bind_rows(CS.rows %>% bind_cols(P.Perm=P.Perms))
  }

  if ("SA" %in% Comps) {
    require(fixest)
    ## Ensure there is a never-treated group:
    Max_Pd <- max(Data$Period)
    Max_ST <- max(Data$Start)
    if (Max_ST <= Max_Pd) {
      Data.SA <- Data %>% filter(Period < Max_ST) %>%
        mutate(Start=if_else(Start==Max_ST, 0, Start))
    }
    SA <- feols(Summ~sunab(cohort=Start,
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
      if ("W_SA.W_ATT" %in% Results$Estimator) {
        SA.row <- SA.row %>% dplyr::mutate(P.Perm=Results %>%
                                             dplyr::filter(Estimator=="W_SA.W_ATT",Outcome=="Probability") %>%
                                             pull(P))
      } else {
        SA.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) feols(Summ~sunab(cohort=Start,
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
    require(DIDmultiplegt)
    CH <- did_multiplegt(
      mode="dyn",
      df=Data,
      outcome="Summ",
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
      if ("W_CH.W_M" %in% Results$Estimator) {
        CH.row <- CH.row %>% dplyr::mutate(P.Perm=Results %>%
                                             dplyr::filter(Estimator=="W_SA.W_ATT",Outcome=="Probability") %>%
                                             pull(P))
      } else {
        CH.Perms <- simplify2array(lapply(Data.Perms.Short,
                                          FUN=function(x) (did_multiplegt(
                                            mode="dyn",
                                            df=Data,
                                            outcome="Summ",
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
    MEM <- lmer(Outcome~Interv+PeriodF+(1|ClusterF), data=Data.LME)
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
    CPI <- lmer(Outcome~Interv+PeriodF+(1|ClusterF)+(1|CPI), data=Data.LME,
                control=lmerControl(check.rankX="silent.drop.cols",
                                    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    CPI.ci <- confint(CPI, parm="Interv", level=CI.Perc, quiet=TRUE)
    CPI.row <- tibble_row(Method="CPI",
                          Estimate=summary(CPI)$coefficients["Interv","Estimate"],
                          SE=summary(CPI)$coefficients["Interv","Std. Error"],
                          CIL=CPI.ci["Interv",1],
                          CIU=CPI.ci["Interv",2]) %>%
      dplyr::mutate(P=2*pnorm(abs(Estimate/SE), lower.tail=FALSE))
    if ("CPI" %in% Comps_PermPs) {
      CPI.Perms <- simplify2array(lapply(Data.Perms,
                                         FUN=function(x) summary(lmer(Outcome~Interv+PeriodF+(1|ClusterF)+(1|CPI), data=x,
                                                                        control=lmerControl(check.rankX="silent.drop.cols",
                                                                                            check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))$coefficients["Interv","Estimate"]))
      CPI.row <- CPI.row %>% dplyr::mutate(P.Perm=mean(abs(CPI.Perms) >= abs(CPI.row$Estimate)))
    }
    Comp_Outs <- Comp_Outs %>% bind_rows(CPI.row)
  }

  if ("CPI.T" %in% Comps) {
    CPI.T <- lmer(Outcome~Interv+Interv:PeriodF+PeriodF+(1|ClusterF)+(1|CPI), data=Data.LME,
                  control=lmerControl(check.rankX="silent.drop.cols",
                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    IntVec <- grepl("Interv:",rownames(summary(CPI.T)$coefficients))
    Names <- c("CPI.T.TAvg","CPI.T.Interv:Baseline")
    IVec <- matrix(c((rownames(summary(CPI.T)$coefficients)=="Interv")+IntVec/(sum(IntVec)+1),
                     as.numeric(rownames(summary(CPI.T)$coefficients)=="Interv")),
                   nrow=2, byrow=TRUE)
    for (i in 1:length(IntVec)) {
      if (IntVec[i]) {
        IVec <- rbind(IVec, as.numeric(rownames(summary(CPI.T)$coefficients)=="Interv") + c(rep(0,i-1),1,rep(0,length(IntVec)-i)))
        Names <- c(Names,paste0("CPI.T.",rownames(summary(CPI.T)$coefficients)[i]))
      }
    }

    CPI.T.row <- tibble(Method=Names,
                            Estimate=(IVec %*% summary(CPI.T)$coefficients[,"Estimate"])[,1],
                            SE=diag(sqrt(IVec %*% vcov(CPI.T) %*% t(IVec)))) %>%
      dplyr::mutate(P=2*pnorm(abs(Estimate/SE), lower.tail=FALSE),
                    CIL=Estimate+qnorm((1-CI.Perc)/2)*SE,
                    CIU=Estimate+qnorm((1+CI.Perc)/2)*SE)
    if ("CPI.T" %in% Comps_PermPs) {
      CPI.T.Perms <- do.call("cbind",
                             lapply(Data.Perms,
                                           FUN=function(x) abs(IVec %*% summary(lmer(Outcome~Interv+Interv:PeriodF+PeriodF+(1|ClusterF)+(1|CPI), data=x,
                                                                                    control=lmerControl(check.rankX="silent.drop.cols",
                                                                                                        check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))$coefficients[,"Estimate"]) >= abs(CPI.T.row$Estimate)))
      CPI.T.row <- CPI.T.row %>% dplyr::mutate(P.Perm=apply(CPI.T.Perms, 1, mean))
    }
    Comp_Outs <- Comp_Outs %>% bind_rows(CPI.T.row)
  }

  if ("CPI.D" %in% Comps) {
    CPI.D <- lmer(Outcome~Interv+Interv:DiffF+PeriodF+(1|ClusterF)+(1|CPI), data=Data.LME,
                  control=lmerControl(check.rankX="silent.drop.cols",
                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    IntVec <- grepl("Interv:",rownames(summary(CPI.D)$coefficients))
    Names <- c("CPI.D.DAvg","CPI.D.Interv:Baseline")
    IVec <- matrix(c((rownames(summary(CPI.D)$coefficients)=="Interv")+IntVec/(sum(IntVec)+1),
                     as.numeric(rownames(summary(CPI.D)$coefficients)=="Interv")),
                   nrow=2, byrow=TRUE)
    for (i in 1:length(IntVec)) {
      if (IntVec[i]) {
        IVec <- rbind(IVec, as.numeric(rownames(summary(CPI.D)$coefficients)=="Interv") + c(rep(0,i-1),1,rep(0,length(IntVec)-i)))
        Names <- c(Names,paste0("CPI.D.",rownames(summary(CPI.D)$coefficients)[i]))
      }
    }
    CPI.D.row <- tibble(Method=Names,
                        Estimate=(IVec %*% summary(CPI.D)$coefficients[,"Estimate"])[,1],
                        SE=diag(sqrt(IVec %*% vcov(CPI.D) %*% t(IVec)))) %>%
      dplyr::mutate(P=2*pnorm(abs(Estimate/SE), lower.tail=FALSE),
                    CIL=Estimate+qnorm((1-CI.Perc)/2)*SE,
                    CIU=Estimate+qnorm((1+CI.Perc)/2)*SE)
    if ("CPI.D" %in% Comps_PermPs) {
      CPI.D.Perms <- do.call("cbind",
                             lapply(Data.Perms,
                                    FUN=function(x) abs(IVec %*% summary(lmer(Outcome~Interv+Interv:DiffF+PeriodF+(1|ClusterF)+(1|CPI), data=x,
                                                                              control=lmerControl(check.rankX="silent.drop.cols",
                                                                                                  check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))$coefficients[,"Estimate"]) >= abs(CPI.D.row$Estimate)))
      CPI.D.row <- CPI.D.row %>% dplyr::mutate(P.Perm=apply(CPI.D.Perms, 1, mean))
    }
    Comp_Outs <- Comp_Outs %>% bind_rows(CPI.D.row)
  }

  if ("CPI.DT" %in% Comps) {
    CPI.DT <- lmer(Outcome~Interv+Interv:DiffF:PeriodF+PeriodF+(1|ClusterF)+(1|CPI), data=Data.LME,
                   control=lmerControl(check.rankX="silent.drop.cols",
                                       check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    RN <- rownames(summary(CPI.DT)$coefficients)

    IntVec <- grepl("Interv:",RN)
    BaseVec <- as.numeric(RN=="Interv")
    IVec <- matrix(BaseVec+IntVec/(sum(IntVec)+1),
                   nrow=1, byrow=TRUE)
    Names <- c("CPI.DT.DTAvg")

    Diffs <- Data.LME %>% dplyr::filter(Interv==1) %>%
      dplyr::group_by(DiffF) %>%
      dplyr::summarize(NP = length(unique(PeriodF)))
    for (a in Diffs$DiffF) {
      if (sum(grepl(paste0("Interv:DiffF",a),RN)) > 0) {
        IVec <- rbind(IVec, BaseVec+grepl(paste0("Interv:DiffF",a),RN)/Diffs$NP[Diffs$DiffF==a])
        Names <- c(Names,paste0("CPI.DT.Diff",a))
      }
    }

    Pds <- Data.LME %>% dplyr::filter(Interv==1) %>%
      dplyr::group_by(PeriodF) %>%
      dplyr::summarize(ND = length(unique(DiffF)))
    for (j in Pds$PeriodF) {
      if (sum(grepl("Interv:",RN)*grepl(paste0(":PeriodF",j),RN)) > 0) {
        IVec <- rbind(IVec, BaseVec+grepl("Interv:",RN)*grepl(paste0(":PeriodF",j),RN)/Pds$ND[Pds$PeriodF==j])
        Names <- c(Names,paste0("CPI.DT.Pd",j))
      }
    }

    IVec <- rbind(IVec,
                  apply(IVec[grepl("CPI.DT.Diff",Names),], MARGIN=2, mean),
                  apply(IVec[grepl("CPI.DT.Pd",Names),], MARGIN=2, mean),
                  apply(IVec[grepl("CPI.DT.Pd",Names) & (!grepl(paste0("CPI.DT.Pd",Pds$PeriodF[dim(Pds)[1]]),Names)),],
                        MARGIN=2, mean),
                  as.numeric(RN=="Interv") + as.numeric(grepl("Interv:",RN) & !(grepl(paste0("PeriodF",Pds$PeriodF[dim(Pds)[1]]),RN)))/sum(Pds$ND[Pds$PeriodF != Pds$PeriodF[dim(Pds)[1]]]))
    Names <- c(Names, "CPI.DT.DAvg", "CPI.DT.TAvg", "CPI.DT.TAvgExLast", "CPI.DT.DTAvgExLastT")

    CPI.DT.row <- tibble(Method=Names,
                             Estimate=(IVec %*% summary(CPI.DT)$coefficients[,"Estimate"])[,1],
                             SE=diag(sqrt(IVec %*% vcov(CPI.DT) %*% t(IVec)))) %>%
      dplyr::mutate(P=2*pnorm(abs(Estimate/SE), lower.tail=FALSE),
                    CIL=Estimate+qnorm((1-CI.Perc)/2)*SE,
                    CIU=Estimate+qnorm((1+CI.Perc)/2)*SE)
    if ("CPI.DT" %in% Comps_PermPs) {
      CPI.DT.Perms <- do.call("cbind",
                             lapply(Data.Perms,
                                    FUN=function(x) abs(IVec %*% summary(lmer(Outcome~Interv+Interv:DiffF:PeriodF+PeriodF+(1|ClusterF)+(1|CPI), data=x,
                                                                              control=lmerControl(check.rankX="silent.drop.cols",
                                                                                                  check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))))$coefficients[,"Estimate"]) >= abs(CPI.DT.row$Estimate)))
      CPI.DT.row <- CPI.DT.row %>% dplyr::mutate(P.Perm=apply(CPI.DT.Perms, 1, mean))
    }
    Comp_Outs <- Comp_Outs %>% bind_rows(CPI.DT.row)
  }

  if (("CLWP" %in% Comps) | ("CLWPA" %in% Comps)) {
    if ("MEM" %in% Comps) {
      start_theta <- MEM.row$Estimate
      start_sigma <- sqrt(as.data.frame(VarCorr(MEM))[2,"vcov"]*2)
    } else if ("CPI" %in% Comps) {
      start_theta <- CPI.row$Estimate
      start_sigma <- sqrt(as.data.frame(VarCorr(CPI))[2,"vcov"]*2)
    } else {
      MEM <- lmer(Outcome~Interv+PeriodF+(1|ClusterF), data=Data.LME)
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

  return(Comp_Outs %>% dplyr::mutate(CIW=CIU-CIL))
}
