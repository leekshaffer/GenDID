#######################################
###### File: ToyExample.R #############
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/18 #############
###### Updated 2024/05/07 #############
#######################################

source("A_Const.R")
source("D_F_Const.R")
source("Rank_Analysis.R")
source("Sigmas.R")
source("Solver.R")
source("Var_Min.R")

## Set setting
N <- 2
J <- 3
Clusters <- c("A","B")
StartPeriods <- c(2,3)

A_mat <- gen_A(N,J)

### Assumption (5): Homogeneity ###
v.5 <- 1

DFT_list.5 <- gen_DFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=5)
Solve.5 <- solve_WA(DFT_obj=DFT_list.5,
                    A_mat=A_mat,
                    v=v.5,
                    rank_obj=rank_an(DFT_list.5,v.5),
                    DID_full=TRUE)
MVar.5 <- min_var(solve_obj=Solve.5,
                  A_mat=A_mat,
                  Sigma=diag(1, nrow=ncol(A_mat)))
## Alternate variance options:
MVar.5.CS <- min_var(solve_obj=Solve.5,
                     A_mat=A_mat,
                     Sigma=create_Sigma_CS(rho=0.1, N=2, J=3))
MVar.5.AR <- min_var(solve_obj=Solve.5,
                     A_mat=A_mat,
                     Sigma=create_Sigma_AR1(rho=0.17, N=2, J=3))


### Assumption (4): Calendar-Time Heterogeneity ###
v.4 <- c(1/2,1/2)
v.4.alt <- c(1,0)
v.4.alt.2 <- c(0,1)

DFT_list.4 <- gen_DFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=4)
Solve.4 <- solve_WA(DFT_obj=DFT_list.4,
                    A_mat=A_mat,
                    v=cbind(v.4,v.4.alt,v.4.alt.2),
                    rank_obj=rank_an(DFT_list.4,cbind(v.4,v.4.alt,v.4.alt.2)),
                    DID_full=TRUE)
MVar.4 <- min_var(solve_obj=Solve.4,
                  A_mat=A_mat,
                  Sigma=diag(1, nrow=ncol(A_mat)))

### Assumption (3): Exposure-Time Heterogeneity ###
v.3 <- c(1/2,1/2)
v.3.alt <- c(1,0)

DFT_list.3 <- gen_DFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=3)
Solve.3 <- solve_WA(DFT_obj=DFT_list.3,
                    A_mat=A_mat,
                    v=cbind(v.3,v.3.alt),
                    rank_obj=rank_an(DFT_list.3,cbind(v.3,v.3.alt)),
                    DID_full=TRUE)
MVar.3 <- min_var(solve_obj=Solve.3,
                  A_mat=A_mat,
                  Sigma=diag(1, nrow=ncol(A_mat)))

### Assumption (2): Exposure- and Calendar-Time Heterogeneity ###
v.2.1 <- c(1,0,0)
v.2.2 <- c(0,1,0)
v.2.3 <- c(0,0,1)

DFT_list.2 <- gen_DFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=2)
Solve.2 <- solve_WA(DFT_obj=DFT_list.2,
                    A_mat=A_mat,
                    v=cbind(v.2.1,v.2.2,v.2.3),
                    rank_obj=rank_an(DFT_list.2,cbind(v.2.1,v.2.2,v.2.3)),
                    DID_full=TRUE)
MVar.2 <- min_var(solve_obj=Solve.2,
                  A_mat=A_mat,
                  Sigma=diag(1, nrow=ncol(A_mat)))


######
### Second setting: N=4,J=3
N <- 4
J <- 3
Clusters <- c("A","B","C","D")
StartPeriods <- c(2,2,3,3)

A_mat <- gen_A(N,J)

### Assumption (5): Homogeneity ###
v.5 <- 1

DFT_list.5 <- gen_DFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=5)
Solve.5 <- solve_WA(DFT_obj=DFT_list.5,
                    A_mat=A_mat,
                    v=v.5,
                    rank_obj=rank_an(DFT_list.5,v.5),
                    DID_full=TRUE)
MVar.5 <- min_var(solve_obj=Solve.5,
                  A_mat=A_mat,
                  Sigma=diag(1, nrow=ncol(A_mat)))

### Assumption (3): Exposure-Time Heterogeneity ###
v.3 <- c(1/2,1/2)
v.3.alt <- c(1,0)

DFT_list.3 <- gen_DFT(Clusters, StartPeriods,
                      OrderedPds=1:3, Assumption=3)
Solve.3 <- solve_WA(DFT_obj=DFT_list.3,
                    A_mat=A_mat,
                    v=matrix(c(.5,.5,1,0), ncol=2, nrow=2),
                    rank_obj=rank_an(DFT_list.3,matrix(c(.5,.5,1,0), ncol=2, nrow=2)),
                    DID_full=TRUE)
MVar.3 <-  min_var(solve_obj=Solve.3,
                   A_mat=A_mat,
                   Sigma=diag(1, nrow=ncol(A_mat)))
