#######################################
###### File: SimpleExample.R ##########
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/18 #############
###### Updated 2024/04/18 #############
#######################################

source("A_Const.R")
source("D_F_Const.R")
source("Rank_Analysis.R")
source("Solver.R")

## Set setting
N <- 2
J <- 3
Clusters <- c("A","B")
StartPeriods <- c(2,3)

A_mat <- gen_A(N,J)

### Assumption (5): Homogeneity ###
v.5 <- 1

DFT_list.5 <- gen_DFT(Clusters, StartPeriods, J, Assumption=5)
Rank_Res.5 <- rank_an(DFT_list.5, v.5)
Solve.5 <- solve_WA(DFT_obj=DFT_list.5,
                    A_mat=A_mat,
                    v=v.5,
                    DID_full=TRUE)

### Assumption (4): Calendar-Time Heterogeneity ###
v.4 <- c(1/2,1/2)

DFT_list.4 <- gen_DFT(Clusters, StartPeriods, J, Assumption=4)
Rank_Res.4 <- rank_an(DFT_list.4, v.4)

v.4.alt <- c(1,0)
Rank_Res.4.alt <- rank_an(DFT_list.4, v.4.alt)

v.4.alt.2 <- c(0,1)
Rank_Res.4.alt.2 <- rank_an(DFT_list.4, v.4.alt.2)

Solve.4 <- solve_WA(DFT_obj=DFT_list.4,
                    A_mat=A_mat,
                    v=cbind(v.4,v.4.alt,v.4.alt.2),
                    DID_full=TRUE)

### Assumption (3): Exposure-Time Heterogeneity ###
v.3 <- c(1/2,1/2)

DFT_list.3 <- gen_DFT(Clusters, StartPeriods, J, Assumption=3)
Rank_Res.3 <- rank_an(DFT_list.3, v.3)

v.3.alt <- c(1,0)
Rank_Res.3.alt <- rank_an(DFT_list.3,v.3.alt)

Solve.3 <- solve_WA(DFT_obj=DFT_list.3,
                    A_mat=A_mat,
                    v=cbind(v.3,v.3.alt),
                    DID_full=TRUE)
