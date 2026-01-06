#######################################
###### File: ToyExample.R #############
###### Lee Kennedy-Shaffer ############
#######################################

library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

source("R/Mat_Const.R")
source("R/Rank_Analysis.R")
source("R/Solver.R")
source("R/Var_Min.R")
source("R/Sigmas.R")

## Set setting
N <- 2
J <- 3
Clusters <- c("A","B")
StartPeriods <- c(2,3)

### Assumption (5): Homogeneity ###
v.5 <- 1

### Intermediate Results:
ADFT_list.5 <- gen_ADFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=5)
ADFT_list.5$A_mat

Solve.5 <- Solve_Assumption(StartTimes=tibble(Cluster=Clusters,
                                              StartPd=StartPeriods),
                            OrderedPds=1:3,
                            Assumption=5,
                            v.Mat=v.5)
Solve.5$DID.weights
Solve.5$Obs.weights

### Main Results:
MVar.5 <- min_var(SolveOut=Solve.5,
                  Sigma=diag(1, nrow=N*J))

MVar.5$MV

## Alternate variance options:
MVar.5.CS <- min_var(SolveOut=Solve.5,
                     Sigma=create_Sigma_CS(rho=0.1, N=2, J=3))
MVar.5.AR <- min_var(SolveOut=Solve.5,
                     Sigma=create_Sigma_AR1(rho=0.17, N=2, J=3))

MVar.5.CS$MV
MVar.5.AR$MV


### Assumption (4): Calendar-Time Heterogeneity ###
v.4 <- c(1/2,1/2)
v.4.alt <- c(1,0)
v.4.alt.2 <- c(0,1)

ADFT_list.4 <- gen_ADFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=4)
Solve.4 <- Solve_Assumption(StartTimes=tibble(Cluster=Clusters,
                                              StartPd=StartPeriods),
                            OrderedPds=1:3,
                            Assumption=4,
                            v.Mat=cbind(v.4,v.4.alt,v.4.alt.2))
MVar.4 <- min_var(SolveOut=Solve.4,
                  Sigma=diag(1, nrow=N*J))

MVar.4$MV

### Assumption (3): Exposure-Time Heterogeneity ###
v.3 <- c(1/2,1/2)
v.3.alt <- c(1,0)

ADFT_list.3 <- gen_ADFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=3)
Solve.3 <- Solve_Assumption(StartTimes=tibble(Cluster=Clusters,
                                              StartPd=StartPeriods),
                            OrderedPds=1:3,
                            Assumption=3,
                            v.Mat=cbind(v.3,v.3.alt))

MVar.3 <- min_var(SolveOut=Solve.3,
                  Sigma=diag(1, nrow=N*J))

MVar.3$MV

### Assumption (2): Exposure- and Calendar-Time Heterogeneity ###
v.2.1 <- c(1,0,0)
v.2.2 <- c(0,1,0)
v.2.3 <- c(0,0,1)

ADFT_list.2 <- gen_ADFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=2)
Solve.2 <- Solve_Assumption(StartTimes=tibble(Cluster=Clusters,
                                              StartPd=StartPeriods),
                            OrderedPds=1:3,
                            Assumption=2,
                            v.Mat=cbind(v.2.1,v.2.2,v.2.3))
MVar.2 <- min_var(SolveOut=Solve.2,
                  Sigma=diag(1, nrow=N*J))

MVar.2$MV


######
### Second setting: N=4,J=3
N <- 4
J <- 3
Clusters <- c("A","B","C","D")
StartPeriods <- c(2,2,3,3)

### Assumption (5): Homogeneity ###
v.5 <- 1

ADFT_list.5 <- gen_ADFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=5)
Solve.5 <- Solve_Assumption(StartTimes=tibble(Cluster=Clusters,
                                              StartPd=StartPeriods),
                            OrderedPds=1:J,
                            Assumption=5,
                            v.Mat=v.5)
MVar.5 <- min_var(SolveOut=Solve.5,
                  Sigma=diag(1, nrow=N*J))

MVar.5$MV

### Assumption (3): Exposure-Time Heterogeneity ###
v.3 <- c(1/2,1/2)
v.3.alt <- c(1,0)

ADFT_list.3 <- gen_ADFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=3)
Solve.3 <- Solve_Assumption(StartTimes=tibble(Cluster=Clusters,
                                              StartPd=StartPeriods),
                            OrderedPds=1:J,
                            Assumption=3,
                            v.Mat=matrix(c(.5,.5,1,0), ncol=2, nrow=2))
MVar.3 <-  min_var(SolveOut=Solve.3,
                   Sigma=diag(1, nrow=N*J))

MVar.3$MV
