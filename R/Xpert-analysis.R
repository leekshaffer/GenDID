#######################################
###### File: Xpert-analysis.R #########
#######################################

## Note: if you are getting errors from some external packages,
## you may wish to leave out "SA" and "CH" from the
## Comps and Comps_PermPs arguments of the Ext_Comps function
## (or omit the Ext_Comps call altogether).

library(tidyverse)
library(lme4)

source("R/Analysis.R")
source("R/Sigmas.R")
source("R/Ext_Comps.R")

## Set seed for reproducibility:
set.seed(20260107)

## Read in (simulated) data from data folder:
load("data/Xpert-Public-Data.Rda")

## View Theta Schematics for Desired Settings
STs <- xpert.dat %>% dplyr::group_by(Cluster) %>%
  dplyr::summarize(StartPd=min(Period[Interv==1]))

### Assumption 5:
A5_Schem <- gen_Theta(js_Obj=gen_js(Clusters=STs$Cluster,
                        StartPeriods=STs$StartPd,
                        OrderedPds=1:max(xpert.dat$Period)),
          Assumption=5)$Schematic
A5_Schem

### Assumption 4:
gen_Theta(js_Obj=gen_js(Clusters=STs$Cluster,
                        StartPeriods=STs$StartPd,
                        OrderedPds=1:max(xpert.dat$Period)),
          Assumption=4)$Schematic

### Assumption 3:
gen_Theta(js_Obj=gen_js(Clusters=STs$Cluster,
                        StartPeriods=STs$StartPd,
                        OrderedPds=1:max(xpert.dat$Period)),
          Assumption=3)$Schematic

### Assumption 2:
gen_Theta(js_Obj=gen_js(Clusters=STs$Cluster,
                        StartPeriods=STs$StartPd,
                        OrderedPds=1:max(xpert.dat$Period)),
          Assumption=2)$Schematic

## Overall Schematic Figure:
Schem_Dat <- as_tibble(A5_Schem) %>%
  dplyr::mutate(Cluster=row_number()) %>%
  tidyr::pivot_longer(cols=-Cluster, names_to="Period") %>%
  dplyr::mutate(Treatment=factor(value,
                                 levels=c(0,1),
                                 labels=c("Control","Intervention")),
                Period=as.numeric(Period))
ggsave(filename="figs/Xpert_Schematic.eps",
       plot=ggplot(data=Schem_Dat,
                   mapping=aes(x=Period, y=Cluster, fill=Treatment)) +
         geom_tile(color="grey80", lty=1) + theme_bw() +
         coord_cartesian(xlim=c(0.5,dim(A5_Schem)[2]+0.5), ylim=c(dim(A5_Schem)[1]+0.5,0.5),
                         clip="off", expand=FALSE) +
         scale_y_reverse(breaks=1:(dim(A5_Schem)[1]), minor_breaks=NULL) +
         scale_x_continuous(breaks=1:(dim(A5_Schem)[2]), minor_breaks=NULL) +
         scale_fill_manual(values=c("white","grey20")) +
         labs(x="Period (Month of Study)", y="Cluster", fill="Treatment",
              title="Trial Schematic for Example SWT"),
       width=6, height=4, units="in")

## Prep for Full Analyses:

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

### Variances To Check
#### Using rho = 0.003 from Thompson et al. 2018
#### Alternative with rho = 0.333 to test robustness to input ICC
#### Alternative with AR1 structure (rho = 0.012 gives average ICC within a cluster ~0.003):
N <- length(unique(xpert.dat$Cluster))
J <- length(unique(xpert.dat$Period))

SigmaList <- list(Ind=create_Sigma_Ind(N=N, J=J),
                  CS_0_003=create_Sigma_CS(rho=0.003, N=N, J=J),
                  CS_0_333=create_Sigma_CS(rho=0.333, N=N, J=J),
                  AR1_0_012=create_Sigma_AR1(rho=0.012, N=N, J=J))

### Values To Check for Confidence Intervals
CIs_Check <- list(Zero=tibble(Probability=0, `Log Odds`=0),
                  NP=tibble(Probability=-0.048, `Log Odds`=log(0.78)),
                  PWP=tibble(Probability=-0.042, `Log Odds`=log(0.85)),
                  NPneg=tibble(Probability=0.048, `Log Odds`=log(1/0.78)))
CIs_Mesh_Vals <- tibble(ValNum=1:701,
                   Name=paste0("SV_",ValNum),
                   Probability=seq(-0.4, 0.3, length.out=701),
                   `Log Odds`=log((Probability+0.5)/(1-Probability-0.5)))
Perms_Use <- 1000

## Full Analysis for Assumption 5:

Analysis.5 <- An_Full(Data=xpert.dat %>% dplyr::rename(Probability=Outcome,
                                                       `Log Odds`=logOdds),
                      OrderedPds=1:max(xpert.dat$Period),
                      Outcome_Cols=c("Probability","Log Odds"),
                      Assumption=5,
                      v.Mat=Vectors.5,
                      SigmaList=SigmaList,
                      Nested_Comps=c("TW","CS","SA","CH","CO","NP"),
                      NumPerms=Perms_Use,
                      CIs_List=CIs_Check,
                      CI.SV.Mesh=CIs_Mesh_Vals,
                      CI.Perc=0.95,
                      Keep_MVO=TRUE)

## Full Analysis for Assumption 4:

Analysis.4 <- An_Full(Data=xpert.dat %>% dplyr::rename(Probability=Outcome,
                                                       `Log Odds`=logOdds),
                      OrderedPds=1:max(xpert.dat$Period),
                      Outcome_Cols=c("Probability","Log Odds"),
                      Assumption=4,
                      v.Mat=Vectors.4,
                      SigmaList=SigmaList,
                      Nested_Comps=c("TW","CS","SA","CH","CO","NP"),
                      NumPerms=Perms_Use,
                      CIs_List=CIs_Check,
                      CI.SV.Mesh=CIs_Mesh_Vals,
                      CI.Perc=0.95,
                      Keep_MVO=TRUE)

## Full Analysis for Assumption 3:

Analysis.3 <- An_Full(Data=xpert.dat %>% dplyr::rename(Probability=Outcome,
                                                       `Log Odds`=logOdds),
                      OrderedPds=1:max(xpert.dat$Period),
                      Outcome_Cols=c("Probability","Log Odds"),
                      Assumption=3,
                      v.Mat=Vectors.3,
                      SigmaList=SigmaList,
                      Nested_Comps=c("TW","CS","SA","CH","CO","NP"),
                      NumPerms=Perms_Use,
                      CIs_List=CIs_Check,
                      CI.SV.Mesh=CIs_Mesh_Vals,
                      CI.Perc=0.95,
                      Keep_MVO=TRUE)

## Full Analysis for Assumption 2:

Analysis.2 <- An_Full(Data=xpert.dat %>% dplyr::rename(Probability=Outcome,
                                                       `Log Odds`=logOdds),
                      OrderedPds=1:max(xpert.dat$Period),
                      Outcome_Cols=c("Probability","Log Odds"),
                      Assumption=2,
                      v.Mat=Vectors.2,
                      SigmaList=SigmaList,
                      Nested_Comps=c("TW","CS","SA","CH","CO","NP"),
                      NumPerms=Perms_Use,
                      CIs_List=CIs_Check,
                      CI.SV.Mesh=CIs_Mesh_Vals,
                      CI.Perc=0.95,
                      Keep_MVO=TRUE)

## External Comparison Methods:

### Data Prep:
xpert.dat.long <- NULL
for (i in 1:(dim(xpert.dat)[1])) {
  add <- xpert.dat[i,,drop=FALSE] %>% dplyr::select(Interv,Period,Cluster) %>%
    dplyr::cross_join(tibble(Outcome=c(rep(1,xpert.dat[i,"Events"]),
                                       rep(0,xpert.dat[i,"Indivs"]-xpert.dat[i,"Events"]))))
  xpert.dat.long <- xpert.dat.long %>% dplyr::bind_rows(add)
}
xpert.dat.long <- xpert.dat.long %>%
  dplyr::left_join(STs %>% rename(Start=StartPd),
                   by=join_by(Cluster))

Comp_Results <- Ext_Comps(Data.Long=xpert.dat.long,
                      SummOutName=NULL,
                      Comps=c("TW", "CS", "SA", "CH",
                              "MEM",
                              "CPI", "CPI.T", "CPI.D", "CPI.DT",
                              "CLWP", "CLWPA"),
                      Comps_PermPs=c("TW", "CS", "SA", "CH",
                                     "MEM",
                                     "CPI", "CPI.T", "CPI.D", "CPI.DT",
                                     "CLWP", "CLWPA"),
                      P.Orders=Analysis.5$Permutation.Orders,
                      NumPerms=Perms_Use,
                      Results=Analysis.5$Results,
                      CI.Perc=0.95) %>%
  bind_rows(Ext_Comps(Data.Long=xpert.dat.long %>%
                        left_join(xpert.dat %>% dplyr::select(Cluster,Period,logOdds) %>%
                                    dplyr::rename(Summ=logOdds),
                                  by=join_by(Cluster,Period)),
                      SummOutName="Summ",
                      Comps=c("CLWP","CLWPA"),
                      Comps_PermPs=c("CLWP","CLWPA"),
                      P.Orders=Analysis.5$Permutation.Orders,
                      NumPerms=Perms_Use,
                      Results=NULL,
                      CI.Perc=0.95) %>%
              dplyr::mutate(Method=c("CLWP Log Odds","CLWPA Log Odds")))

## Save Results:

save(list=c("Analysis.5","Analysis.4","Analysis.3","Analysis.2",
            "Comp_Results"),
     file="int/Xpert-Public-Output.Rda")


### Observation Weight Heatmaps:

### To create various heat maps, add rows with
### different values of i (Assumption Setting),
### j (Variance setting), and Estimators (estimator)

Map_Settings <- tibble(i=c(5,4,3,2,3,2,rep(4,6),rep(2,6)),
                       j=rep("CS_0_003",18),
                       Estimator=c("Single","AvgEx8","Avg","AvgExT8",
                                    "D.1","D.1",
                                    "T.2","T.3","T.4","T.5","T.6","T.7",
                                    "T.2","T.3","T.4","T.5","T.6","T.7"),
                       Est_labs=c("Overall, Assumption S5","Avg., Assumption S4",
                                  "Avg., Assumption S3", "ATT, Assumption S2",
                                  "First-Period, Assumption S3", "First-Period, Assumption S2",
                                  "Pd. 2, Assumption S4", "Pd. 3, Assumption S4",
                                  "Pd. 4, Assumption S4", "Pd. 5, Assumption S4",
                                  "Pd. 6, Assumption S4", "Pd. 7, Assumption S4",
                                  "Pd. 2, Assumption S2", "Pd. 3, Assumption S2",
                                  "Pd. 4, Assumption S2", "Pd. 5, Assumption S2",
                                  "Pd. 6, Assumption S2", "Pd. 7, Assumption S2"))
for (row in 1:(dim(Map_Settings)[1])) {
  Weights <- ((get(paste0("Analysis.",Map_Settings$i[row]))[["MVO"]])[[Map_Settings$j[row]]])$MV$Obs.weights
  if (is.null(colnames(Weights))) {
    colnames(Weights) <- as.character(1:(dim(Weights)[2]))
  }
  Obs.weight.dat <- tibble(x=rep(1:J, times=N), y=rep(1:N, each=J),
                           Value=Weights[,Map_Settings$Estimator[row]])
  ggsave(filename=paste0("figs/Xpert-Weights_Heatmap_",Map_Settings[row,"i"],"_",
                         Map_Settings$j[row],"_",
                         Map_Settings$Estimator[row],".eps"),
         plot=ggplot(data=Obs.weight.dat, mapping=aes(x=x, y=y, fill=Value)) +
           geom_tile(color="grey80", lty=1) + theme_bw() +
           geom_text(aes(label=format(round(Value, digits=3), nsmall=3))) +
           coord_cartesian(xlim=c(0.5,J+0.5), ylim=c(N+0.5,0.5), clip="off", expand=FALSE) +
           scale_y_reverse(breaks=1:N, minor_breaks=NULL) +
           scale_x_continuous(breaks=1:J, minor_breaks=NULL) +
           scale_fill_gradient2(low="#7938C6",high="#DE6D07") +
           labs(x="Period (Month of Study)", y="Cluster", fill="Weight",
                title=paste0("Observation Weights: ",
                             Map_Settings$Est_labs[row])),
         width=6, height=4, units="in")
}

## MS Results:

Tbl2 <- Analysis.5$Results %>%
  dplyr::filter(Estimator=="CS_0_003_Single") %>%
  dplyr::mutate(Assumption="S5") %>%
  dplyr::bind_rows(Analysis.4$Results %>%
              dplyr::filter(Estimator=="CS_0_003_AvgEx8") %>%
              dplyr::mutate(Assumption="S4")) %>%
  dplyr::bind_rows(Analysis.3$Results %>%
              dplyr::filter(Estimator=="CS_0_003_Avg") %>%
              dplyr::mutate(Assumption="S3")) %>%
  dplyr::bind_rows(Analysis.2$Results %>%
              dplyr::filter(Estimator=="CS_0_003_AvgExT8") %>%
              dplyr::mutate(Assumption="S2")) %>%
  dplyr::select(Assumption,Estimator,Outcome,Estimate,P,CIL,CIU) %>%
  dplyr::mutate(Outcome=if_else(Outcome=="Probability",
                                "RD",
                                "OR")) %>%
  tidyr::pivot_wider(id_cols=c(Assumption,Estimator),
                     names_from=Outcome,
                     values_from=c(Estimate,P,CIL,CIU)) %>%
  dplyr::mutate(Estimate_RD=100*Estimate_RD,
                CIL_RD=100*CIL_RD,
                CIU_RD=100*CIU_RD,
                Estimate_OR=exp(Estimate_OR),
                CIL_OR=exp(CIL_OR),
                CIU_OR=exp(CIU_OR)) %>%
  dplyr::select(Assumption,Estimator,ends_with("_OR"),ends_with("_RD"))

Vars <- c(Analysis.5$MVO$CS_0_003$MV$Variance[1,"Single"],
          Analysis.4$MVO$CS_0_003$MV$Variance[1,"AvgEx8"],
          Analysis.3$MVO$CS_0_003$MV$Variance[1,"Avg"],
          Analysis.2$MVO$CS_0_003$MV$Variance[1,"AvgExT8"])
Rel_Effs <- c(Vars[2]/Vars[1], Vars[3]/Vars[1], Vars[4]/Vars[1])

Ts <- paste0("CS_0_003_T.",2:7)
Tbl3 <- Analysis.4$Results %>%
  dplyr::filter(Estimator %in% Ts) %>%
  dplyr::mutate(Assumption="S4") %>%
  bind_rows(Analysis.2$Results %>%
              dplyr::filter(Estimator %in% Ts) %>%
              dplyr::mutate(Assumption="S2")) %>%
  dplyr::select(Estimator,Assumption,Outcome,Estimate,P,CIL,CIU) %>%
  dplyr::mutate(Outcome=if_else(Outcome=="Probability",
                                "RD",
                                "OR")) %>%
  tidyr::pivot_wider(id_cols=c(Estimator,Assumption),
                     names_from=c(Outcome),
                     values_from=c(Estimate,P,CIL,CIU)) %>%
  dplyr::mutate(Estimate_RD=100*Estimate_RD,
                CIL_RD=100*CIL_RD,
                CIU_RD=100*CIU_RD,
                Estimate_OR=exp(Estimate_OR),
                CIL_OR=exp(CIL_OR),
                CIU_OR=exp(CIU_OR),
                Month=as.numeric(sub("CS_0_003_T.","",Estimator))) %>%
  tidyr::pivot_wider(id_cols=Month,
                     names_from=c(Assumption),
                     values_from=ends_with(c("_RD","_OR"))) %>%
  dplyr::select(Month,ends_with("OR_S4"),ends_with("OR_S2"),
                ends_with("RD_S4"),ends_with("RD_S2"))

TblWA34 <- Analysis.5$Results %>%
  dplyr::filter(Estimator %in% paste0(c("Ind_","CS_0_003_","AR1_0_012_"),"Single")) %>%
  dplyr::mutate(Assumption="S5") %>%
  dplyr::bind_rows(Analysis.4$Results %>%
                     dplyr::filter(Estimator %in% paste0(c("Ind_","CS_0_003_","AR1_0_012_"),"AvgEx8")) %>%
                     dplyr::mutate(Assumption="S4")) %>%
  dplyr::bind_rows(Analysis.3$Results %>%
                     dplyr::filter(Estimator %in% paste0(c("Ind_","CS_0_003_","AR1_0_012_"),"Avg")) %>%
                     dplyr::mutate(Assumption="S3")) %>%
  dplyr::bind_rows(Analysis.2$Results %>%
                     dplyr::filter(Estimator %in% paste0(c("Ind_","CS_0_003_","AR1_0_012_"),"AvgExT8")) %>%
                     dplyr::mutate(Assumption="S2")) %>%
  dplyr::mutate(Variance=sub("_.*", "", Estimator)) %>%
  dplyr::select(Assumption,Variance,Outcome,Estimate,P,CIL,CIU) %>%
  dplyr::mutate(Outcome=if_else(Outcome=="Probability",
                                "RD",
                                "OR")) %>%
  tidyr::pivot_wider(id_cols=c(Assumption,Variance),
                     names_from=Outcome,
                     values_from=c(Estimate,P,CIL,CIU)) %>%
  dplyr::mutate(Estimate_RD=100*Estimate_RD,
                CIL_RD=100*CIL_RD,
                CIU_RD=100*CIU_RD,
                Estimate_OR=exp(Estimate_OR),
                CIL_OR=exp(CIL_OR),
                CIU_OR=exp(CIU_OR)) %>%
  pivot_wider(id_cols=Assumption,
              names_from=Variance,
              values_from=ends_with(c("_RD","_OR"))) %>%
  dplyr::select(Assumption,ends_with("_Ind"),ends_with("_CS"),ends_with("_AR1"))

TblWA3 <- TblWA34 %>%
  dplyr::select(Assumption,contains("_OR_"))
colnames(TblWA3) <- sub("_OR_","_",colnames(TblWA3))

TblWA4 <- TblWA34 %>%
  dplyr::select(Assumption,contains("_RD_"))
colnames(TblWA4) <- sub("_RD_","_",colnames(TblWA4))

TblWA5C1 <- Analysis.5$Results %>%
  dplyr::filter(Estimator %in% c("W_TW","W_CS.W_simple","W_CS.W_dynamic",
                                 "W_CS.W_group","W_CS.W_calendar",
                                 "W_SA.W_ATT","W_CH.W_M",
                                 "W_CO.W_CO1","W_CO.W_CO2",
                                 "W_CO.W_CO3","W_NP_Eq",
                                 "W_NP_ATT","W_NP_IV"),
                Outcome=="Probability") %>%
  dplyr::select(Estimator,Estimate,P,CIL,CIU) %>%
  bind_rows(Comp_Results %>%
              dplyr::filter(Method %in% c("CPI","CPI.DT","CLWP","CLWPA")) %>%
              dplyr::select(Method,Estimate,P.Perm,CIL,CIU) %>%
              dplyr::rename(Estimator=Method,
                            P=P.Perm))%>%
  dplyr::mutate(Estimate=100*Estimate,
                CIL=100*CIL,
                CIU=100*CIU)

TblWA5C2 <- Analysis.5$Results %>%
  dplyr::filter(Estimator=="CS_0_003_Single",
                Outcome=="Probability") %>%
  bind_rows(Analysis.2$Results %>%
              dplyr::filter(Estimator=="CS_0_003_AvgExT8",
                            Outcome=="Probability")) %>%
  bind_rows(Analysis.3$Results %>%
              dplyr::filter(Estimator=="CS_0_003_Avg",
                            Outcome=="Probability")) %>%
  bind_rows(Analysis.2$Results %>%
              dplyr::filter(Estimator=="CS_0_003_Group",
                            Outcome=="Probability")) %>%
  bind_rows(Analysis.4$Results %>%
              dplyr::filter(Estimator=="CS_0_003_AvgEx8",
                            Outcome=="Probability")) %>%
  bind_rows(Analysis.2$Results %>%
              dplyr::filter(Estimator=="CS_0_003_AvgExT8",
                            Outcome=="Probability")) %>%
  bind_rows(Analysis.2$Results %>%
              dplyr::filter(Estimator=="CS_0_003_D.1",
                            Outcome=="Probability")) %>%
  bind_rows(Analysis.2$Results %>%
              dplyr::filter(Estimator=="CS_0_003_D.1",
                            Outcome=="Probability")) %>%
  bind_rows(Analysis.3$Results %>%
              dplyr::filter(Estimator=="CS_0_003_D.1",
                            Outcome=="Probability")) %>%
  bind_rows(Analysis.5$Results %>%
              dplyr::filter(Estimator=="CS_0_003_Single",
                            Outcome=="Probability")) %>%
  dplyr::select(Estimator,Estimate,P,CIL,CIU) %>%
  dplyr::mutate(Estimate=100*Estimate,
                CIL=100*CIL,
                CIU=100*CIU) %>%
  dplyr::rename_with(.fn=~paste0("GD_",.x))

Padding <- dim(TblWA5C1)[1] - dim(TblWA5C2)[1]
TblWA5 <- TblWA5C1 %>%
  bind_cols(TblWA5C2 %>%
              bind_rows(tibble(GD_Estimator=rep("-",Padding),
                               GD_Estimate=rep(NA,Padding),
                               GD_P=rep(NA,Padding),
                               GD_CIL=rep(NA,Padding),
                               GD_CIU=rep(NA,Padding))))

save(list=c("Tbl2","Vars","Rel_Effs","Tbl3",
            "TblWA3","TblWA4","TblWA5"),
     file="res/Xpert-Public-Results.Rda")


