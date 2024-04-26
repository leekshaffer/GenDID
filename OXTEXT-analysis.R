#######################################
###### File: OXTEXT-analysis.R ########
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
oxtext <- read_xlsx(path="data/oxtext7.xlsx",
                    sheet="OXTEXT 7 Data")

## Get unique periods, clusters, start times, N, J:
Periods <- unique(oxtext$Study_Month)
J <- length(Periods)

StartTimes <- oxtext %>% dplyr::filter(Treatment_Arm=="Treatment") %>% 
  group_by(Centre) %>% dplyr::summarize(StartPd=min(Study_Month)) %>%
  dplyr::arrange(StartPd)
N <- length(StartTimes$Centre)

## Generate A,D,F,Theta matrices:
Amat <- gen_A(N,J)
DFT_5 <- gen_DFT(Clusters=StartTimes$Centre,
                   StartPeriods=StartTimes$StartPd,
                   J=J,
                   Assumption=5)

## Solve:
system.time(solve_5 <- solve_WA(DFT_obj=DFT_5,
                    A_mat=Amat,
                    v=1,
                    rank_obj=NULL, DID_full=TRUE))
## Note: took ~1hr10min on laptop

## Minimize Variance:
# Rho = 0.03 comes from Nickless paper analysis of OXTEXT-7 data
system.time(mv_5_CS <- min_var(solve_obj=solve_5,
                            A_mat=Amat,
                            Sigma=create_Sigma_CS(rho=0.03,N=N,J=J)))
## Note: A5 took ~15sec on laptop

save(list=c("Amat","DFT_5","solve_5","mv_5_CS"), file="../int-large/oxtext-a5.Rda")

## Prep Data:
Ordered_Data <- oxtext %>% dplyr::left_join(DFT_5$Theta$Full %>% 
                                              dplyr::select(Clusters,Cl.Num) %>% distinct(), 
                                     by=join_by(Centre==Clusters)) %>%
  dplyr::group_by(Centre,Cl.Num,Study_Month) %>%
  dplyr::summarize(AvgY=mean(HONOS_Score, na.rm=TRUE), 
                   VarY=var(HONOS_Score, na.rm=TRUE),
                   Num_Patients=n()) %>%
  arrange(Cl.Num,Study_Month)

## Apply Analysis:
Obs_Y <- matrix(Ordered_Data$AvgY, ncol=1)
Obs_D <- Amat %*% Obs_Y
Estimate <- t(mv_5_CS$Obs.weights) %*% Obs_Y

## Put Weights onto D matrix
D_Vals <- DFT_5$D_aug
D_Vals$Value <- Amat %*% Obs_Y
D_Vals$Weight <- mv_5_CS$DID.weights[,1]
D_Full <- D_Vals %>% arrange(desc(Weight))

save(list=c("Amat","Obs_Y","Estimate","D_Full"), file="int/oxtext-a5-res.Rda")
