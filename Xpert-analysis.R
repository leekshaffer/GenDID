#######################################
###### File: Xpert-analysis.R #########
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/25 #############
###### Updated 2024/04/25 #############
#######################################

require(readxl)
require(tidyverse)

source("A_Const.R")
source("D_F_Const.R")
source("Sigmas.R")
source("Rank_Analysis.R")
source("Var_Min.R")
source("Solver.R")

## Read in data:
load("data/Xpert-data.Rda")

## Get unique periods, clusters, start times, N, J:
Periods <- unique(trajclust2$Period)
J <- length(Periods)

StartTimes <- trajclust2 %>% dplyr::filter(Interv==1) %>%
  group_by(Cluster) %>% dplyr::summarize(StartPd=min(Period)) %>%
  dplyr::arrange(StartPd)
N <- length(StartTimes$Cluster)

## Generate A,D,F,Theta matrices:
Amat <- gen_A(N,J)

###
DFT_4 <- gen_DFT(Clusters=StartTimes$Cluster,
                 StartPeriods=StartTimes$StartPd,
                 J=J,
                 Assumption=4)
v.4 <- c(rep(1/6,6),0)

## Solve:
system.time(rank_4 <- rank_an(DFT_obj=DFT_4,
                              v=v.4))
## Note: A4 took ~0.2sec on laptop
system.time(solve_4 <- solve_WA(DFT_obj=DFT_4,
                                A_mat=Amat,
                                v=v.4,
                                rank_obj=rank_4, DID_full=TRUE))
## Note: A5 took ~3min on laptop, about the same amount of time with DID_full True or False
## Note: A4 took ~1min on desktop, ~3min on laptop

## Minimize Variance:
# Rho = 0.003 comes from Thompson paper analysis of OXTEXT-7 data
system.time(mv_4_CS <- min_var(solve_obj=solve_4,
                               A_mat=Amat,
                               Sigma=create_Sigma_CS(rho=0.003,N=N,J=J)))
## Note: took ~1.5sec on laptop, .5 on desktop
## Longer version (full DID vectors) for A4: took 37sec on desktop

save(list=c("Amat","DFT_4","solve_4","mv_4_CS"), file="../int_large/xpert-a4.Rda")

## Prep Data:
Ordered_Data_4 <- trajclust2 %>% dplyr::left_join(DFT_4$Theta$Full %>% 
                                              dplyr::select(Clusters,Cl.Num) %>% distinct(), 
                                            by=join_by(Cluster==Clusters)) %>%
  arrange(Cl.Num,Period)

## Apply Analysis:
Obs_Y <- matrix(data=c(Ordered_Data_4$Outcome, Ordered_Data_4$logOdds), ncol=2)
colnames(Obs_Y) <- c("Probability","Log Odds")
Estimate <- t(mv_4_CS$Obs.weights) %*% Obs_Y
Estimate <- c(Estimate,exp(Estimate[1,"Log Odds"]))
names(Estimate) = c("Risk Difference","Log Odds Ratio","Odds Ratio")

## Put Weights onto D matrix
D_Vals <- DFT_4$D_aug
D_Vals$RD <- Amat %*% Obs_Y[,"Probability"]
D_Vals$LOR <- Amat %*% Obs_Y[,"Log Odds"]
D_Vals$OR <- exp(D_Vals$LOR)
D_Vals$Weight <- mv_4_CS$DID.weights[,1]
D_Full <- D_Vals %>% arrange(desc(Weight))
save(list=c("Amat","Obs_Y","Estimate","D_Full"), file="int/xpert-a4-res.Rda")
