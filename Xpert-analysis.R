#######################################
###### File: Xpert-analysis.R #########
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/25 #############
###### Updated 2024/08/08 #############
#######################################

require(readxl)
require(tidyverse)

source("A_Const.R")
source("Sigmas.R")
source("Full_Analysis.R")
source("CompEsts.R")

## Read in (simulated) data from data folder:
load("data/Xpert-data-sim.Rda")

## Get unique periods, clusters, start times, N, J:
Periods <- unique(xpert.dat$Period)
OrderedPds <- Periods[order(Periods)]
J <- length(OrderedPds)

StartTimes <- xpert.dat %>% dplyr::filter(Interv==1) %>%
  group_by(Cluster) %>% dplyr::summarize(StartPd=min(Period)) %>%
  dplyr::arrange(StartPd,Cluster)
N <- length(StartTimes$Cluster)

## Prep Outcome Data in appropriate order:
Ord_Data <- xpert.dat %>% left_join(StartTimes, by="Cluster") %>%
  dplyr::arrange(StartPd,Cluster,Period)
Obs_Y <- matrix(data=c(Ord_Data$Outcome, Ord_Data$logOdds), ncol=2)
colnames(Obs_Y) <- c("Probability","Log Odds")

## Generate A matrix:
Amat <- gen_A(N,J)

## Run Solver for different assumption settings:
SO2 <- Solve_Assumption(Amat,StartTimes,OrderedPds,
                        Assumption=2,
                        v.Mat=cbind(Avg=c(rep(1/28,28)),
                                    AvgEx7=c(rep(1/21,21),rep(0,7)),
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
                                    U.21=c(rep(0,20),1,rep(0,7))),
                        save_loc="../int_large/",
                        save_prefix="xpert-solve-a_")
SO3 <- Solve_Assumption(Amat,StartTimes,OrderedPds,
                        Assumption=3,
                        v.Mat=cbind(Avg=rep(1/7,7),
                             AvgEx7=c(rep(1/6,6),0),
                             D.1=c(1,rep(0,6)),
                             D.2=c(0,1,rep(0,5)),
                             D.3=c(0,0,1,rep(0,4)),
                             D.4=c(0,0,0,1,rep(0,3)),
                             D.5=c(rep(0,4),1,0,0),
                             D.6=c(rep(0,5),1,0),
                             D.7=c(rep(0,6),1),
                             Middle=c(0,1/3,1/3,1/3,0,0,0)),
                        save_loc="../int_large/",
                        save_prefix="xpert-solve-a_")
SO4 <- Solve_Assumption(Amat,StartTimes,OrderedPds,
                        Assumption=4,
                        v.Mat=cbind(AvgEx7=c(rep(1/6,6),0),
                             T.2=c(1,rep(0,6)),
                             T.3=c(0,1,rep(0,5)),
                             T.4=c(0,0,1,rep(0,4)),
                             T.5=c(0,0,0,1,rep(0,3)),
                             T.6=c(rep(0,4),1,0,0),
                             T.7=c(rep(0,5),1,0),
                             T.8=c(rep(0,6),1),
                             Middle=c(0,1/3,1/3,1/3,0,0,0)),
                        save_loc="../int_large/",
                        save_prefix="xpert-solve-a_")
SO5 <- Solve_Assumption(Amat,StartTimes,OrderedPds,
                        Assumption=5,
                        v.Mat=1,
                        save_loc="../int_large/",
                        save_prefix="xpert-solve-a_")

set.seed(413354)

## Run variance minimizer for different settings:
### Independence:
for (i in 2:5) {
  assign(x=paste0("MVOut_",i,"_Ind"),
         value=MV_Assumption(SolveOut=get(paste0("SO",i)),
                             Assumption=i,
                             Sigma=create_Sigma_Ind(N=N,J=J),
                             SigmaName="Ind",
                             Observations=Obs_Y,
                             Permutations=1000,
                             save_loc="int/",
                             save_prefix="xpert-mv-a_"))
}

### Exchangeable (rho = 0.003 from Thompson et al. 2018):
for (i in 2:5) {
  assign(x=paste0("MVOut_",i,"_CS_0_003"),
         value=MV_Assumption(SolveOut=get(paste0("SO",i)),
                             Assumption=i,
                             Sigma=create_Sigma_CS(rho=0.003,N=N,J=J),
                             SigmaName="CS_0_003",
                             Observations=Obs_Y,
                             Permutations=1000,
                             save_loc="int/",
                             save_prefix="xpert-mv-a_"))
}

### AR(1) (rho = 0.012 gives average ICC within a cluster ~0.003):
for (i in 2:5) {
  assign(x=paste0("MVOut_",i,"_AR1_0_012"),
         value=MV_Assumption(SolveOut=get(paste0("SO",i)),
                             Assumption=i,
                             Sigma=create_Sigma_AR1(rho=0.012,N=N,J=J),
                             SigmaName="AR1_0_012",
                             Observations=Obs_Y,
                             Permutations=1000,
                             save_loc="int/",
                             save_prefix="xpert-mv-a_"))
}

## Import Results:
Assns <- 2:5
SigmaNames <- c("Ind","CS_0_003","AR1_0_012")

for (j in SigmaNames) {
  for (i in Assns) {
    load(file=paste0("int/xpert-mv-a_",i,"_",j,".Rda"))
  }
}

## Summarize Results:
for (j in SigmaNames) {
  print(paste0("Variance ",j))
  for (i in Assns) {
    print(paste0("Assumption ",i))
    print((get(paste0("MVOut_",i,"_",j))[["MV"]])[["Variance"]])
  }
}

for (i in Assns) {
  print(paste0("Assumption ",i))
  ProbEsts <- NULL
  OREsts <- NULL
  for (j in SigmaNames) {
    ProbEsts <- cbind(ProbEsts,(get(paste0("MVOut_",i,"_",j))[["Estimates"]])[,"Probability"])
    OREsts <- cbind(OREsts,exp((get(paste0("MVOut_",i,"_",j))[["Estimates"]])[,"Log Odds"]))
  }
  colnames(ProbEsts) <- SigmaNames
  colnames(OREsts) <- SigmaNames
  print(ProbEsts)
  print(OREsts)
}


### Inference:
for (i in Assns) {
  print(paste0("P-Values for Assumption ",i))
  for (j in SigmaNames) {
    print(j)
    print((get(paste0("MVOut_",i,"_",j)))[["P_Values"]])
  }
}

### Observation Weight Heatmaps:
# for (i in Assns) {
#   for (j in SigmaNames) {
## Uncomment preceding two lines to get all heat maps

i <- 5
j <- "CS_0_003"
    Weights <- (get(paste0("MVOut_",i,"_",j))[["MV"]])[["Obs.weights"]]
    for (n in 1:ncol(Weights)) {
      Obs.weight.dat <- tibble(x=rep(1:J, times=N), y=rep(1:N, each=J),
                               Value=Weights[,n])
      ggsave(filename=paste0("figs/Weights_Heatmap_",i,"_",j,"_Col",n,".png"),
             plot=ggplot(data=Obs.weight.dat, mapping=aes(x=x, y=y, fill=Value)) +
               geom_tile() + theme_bw() + 
               coord_cartesian(xlim=c(0.5,8.5), ylim=c(14.5,0.5), clip="off", expand=FALSE) +
               scale_y_reverse(breaks=1:N, minor_breaks=NULL) +
               scale_x_continuous(breaks=1:J, minor_breaks=NULL) +
               scale_fill_gradient2(low="#542788",high="#b35806") +
               labs(x="Period", y="Cluster", fill="Weight",
                    title=paste0("Observation Weights, Assumption: ",i)),
             width=6, height=4, units="in", dpi=600)
    }
#   }
# }
    
    
## Comparisons to other methods:
DFT <- SO5$DFT
Comp_wts <- Comp_Ests_Weights(DFT_obj=DFT, Amat=Amat,
                              estimator=c("CS","SA","CH","CO","NP"))
Comp_ests <- t(as.matrix(Comp_wts$Obs.weights)) %*% Obs_Y
Comp_ests
save(Comp_ests, file="int/Comparison_estimates.Rda")


## Check against existing packages for staggered adoption methods:
### Packages:
require(did) ## For CS
require(fixest) ## For SA
require(DIDmultiplegt) ## For CH

### Data Prep:
xpert.dat.2 <- xpert.dat %>% left_join(StartTimes, by="Cluster") %>%
  mutate(ClusterF=factor(Cluster),
         LeadLag=Period-StartPd)
#### A version without Period 8 to create a never-treated group:
xpert.dat.2.ex8 <- xpert.dat.2  %>% filter(Period != 8) %>%
  mutate(StartPd=if_else(StartPd==8,0,StartPd))

### Callaway and Sant'Anna (2021):
CS_gt <- att_gt(yname="Outcome",
                tname="Period",
                idname="Cluster",
                gname="StartPd",
                data=xpert.dat.2,
                panel=TRUE,
                control_group="notyettreated")
ggdid(CS_gt)
CS_gt
aggte(CS_gt, type="simple")
aggte(CS_gt, type="dynamic")
aggte(CS_gt, type="group")
aggte(CS_gt, type="calendar")

### Sun and Abraham (2021):
SA <- feols(Outcome~sunab(cohort=StartPd,
                          period=Period,
                          ref.c=c(8),
                          ref.p=c(-1,6), att=TRUE) | Cluster + Period,
            data=xpert.dat.2.ex8)
summary(SA)

### de Chaisemartin and d'Haultfoeuille (2020):
CH <- did_multiplegt(df=xpert.dat,
                     Y="Outcome",
                     G="Cluster",
                     T="Period",
                     D="Interv")
CH
