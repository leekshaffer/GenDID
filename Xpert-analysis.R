#######################################
###### File: Xpert-analysis.R #########
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/25 #############
###### Updated 2024/04/26 #############
#######################################

source("A_Const.R")
source("Sigmas.R")
source("Full_Analysis.R")

## Read in data:
load("data/Xpert-data.Rda")

## Get unique periods, clusters, start times, N, J:
Periods <- unique(trajclust2$Period)
J <- length(Periods)

StartTimes <- trajclust2 %>% dplyr::filter(Interv==1) %>%
  group_by(Cluster) %>% dplyr::summarize(StartPd=min(Period)) %>%
  dplyr::arrange(StartPd,Cluster)
N <- length(StartTimes$Cluster)

## Prep Outcome Data in appropriate order:
Ord_Data <- trajclust2 %>% left_join(StartTimes, by="Cluster") %>%
  dplyr::arrange(StartPd,Cluster,Period)
Obs_Y <- matrix(data=c(Ord_Data$Outcome, Ord_Data$logOdds), ncol=2)
colnames(Obs_Y) <- c("Probability","Log Odds")

## Generate A matrix:
Amat <- gen_A(N,J)

## Run Solver for different assumption settings:
SO3 <- Solve_Assumption(Amat,StartTimes,J,
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
                             Middle=c(0,1/3,1/3,1/3,0,0,0)))
SO4 <- Solve_Assumption(Amat,StartTimes,J,
                        Assumption=4,
                        v.Mat=cbind(AvgEx7=c(rep(1/6,6),0),
                             T.1=c(1,rep(0,6)),
                             T.2=c(0,1,rep(0,5)),
                             T.3=c(0,0,1,rep(0,4)),
                             T.4=c(0,0,0,1,rep(0,3)),
                             T.5=c(rep(0,4),1,0,0),
                             T.6=c(rep(0,5),1,0),
                             T.7=c(rep(0,6),1),
                             Middle=c(0,1/3,1/3,1/3,0,0,0)))
SO5 <- Solve_Assumption(Amat,StartTimes,J,
                        Assumption=5,
                        v.Mat=1)
## Timing (desktop) notes: <1min each for A3/A4/A5

## Run variance minimizer for different settings:

### Independence:
for (i in 3:5) {
  assign(x=paste0("MV",i,"_Ind"),
         value=MV_Assumption(SolveOut=get(paste0("SO",i)),
                             Assumption=i,
                             Sigma=create_Sigma_Ind(N=N,J=J),
                             SigmaName="Ind",
                             Observations=Obs_Y))
}

### Exchangeable (rho = 0.003 from Thompson et al. 2018):
for (i in 3:5) {
  assign(x=paste0("MV",i,"_CS_0_003"),
         value=MV_Assumption(SolveOut=get(paste0("SO",i)),
                             Assumption=i,
                             Sigma=create_Sigma_CS(rho=0.003,N=N,J=J),
                             SigmaName="CS_0_003",
                             Observations=Obs_Y))
}
## Timing (desktop) notes: ~3sec each for A3/A4/A5

### AR(1) (rho = 0.012 gives average ICC within a cluster ~0.003):
for (i in 3:5) {
  assign(x=paste0("MV",i,"_AR1_0_012"),
         value=MV_Assumption(SolveOut=get(paste0("SO",i)),
                             Assumption=i,
                             Sigma=create_Sigma_AR1(rho=0.012,N=N,J=J),
                             SigmaName="AR1_0_012",
                             Observations=Obs_Y))
}

