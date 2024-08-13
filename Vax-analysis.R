#######################################
###### File: Vax-analysis.R ###########
###### Lee Kennedy-Shaffer ############
###### Created 2024/08/09 #############
###### Updated 2024/08/13 #############
#######################################

require(tidyverse)

source("A_Const.R")
source("Sigmas.R")
source("Full_Analysis.R")
source("CompEsts.R")

# Fuller data:
## Read in data:
load("data/synth_data_clean.RData")
## Create data set:
### Weeks based on https://ndc.services.cdc.gov/wp-content/uploads/W2021-22.pdf:
MMWR_weeks <- tibble(Period=1:52,
                     First_Day=ymd("2021-01-03")+seq(from=0,by=7,length.out=52),
                     Last_Day=ymd("2021-01-09")+seq(from=0,by=7,length.out=52))
### Lottery dates and references/notes from Fuller et al.
Lott_weeks <- lotteries %>% left_join(MMWR_weeks,
                                      by=join_by(lottery_announced >= First_Day,
                                                 lottery_announced <= Last_Day)) %>%
  dplyr::rename(lottery=lottery_incentive, lott_date=lottery_announced,
                lott_week=Period) %>%
  dplyr::select(state,lottery,lott_week,lott_date)
### Create weekly data set and add lottery weeks:
Vax_wk_2021 <- state_vaccines %>% 
  dplyr::rename(First_18Pop_Pct=Administered_Dose1_Recip_18PlusPop_Pct,
                Complete_18Pop_Pct=Series_Complete_18PlusPop_Pct,
                Cluster=state, Week_End=date) %>%
  left_join(MMWR_weeks %>% dplyr::rename(Week_Start=First_Day), 
            by=join_by(Week_End == Last_Day)) %>%
  dplyr::filter(!is.na(Period)) %>%
  dplyr::select(Cluster,Period,Week_Start,Week_End,
                First_18Pop_Pct,Complete_18Pop_Pct) %>%
  left_join(Lott_weeks,
            by=join_by(Cluster==state)) %>%
  mutate(Interv=if_else(lottery==0,0,if_else(Period >= lott_week, 1, 0)),
         rel_week=if_else(lottery==0,NA,Period-lott_week),
         lott_week=if_else(is.na(lott_week),Inf,lott_week))
### Add previous week
Vax_prev_wk <- Vax_wk_2021 %>% 
  dplyr::select(Cluster,Period,First_18Pop_Pct,Complete_18Pop_Pct) %>%
  dplyr::mutate(Period=Period+1) %>%
  dplyr::rename(Prev_First=First_18Pop_Pct,
                Prev_Complete=Complete_18Pop_Pct)
Vax_wk_2021 <- Vax_wk_2021 %>% 
  left_join(Vax_prev_wk, by=c("Cluster","Period")) %>%
  dplyr::filter(!is.na(Prev_First)) %>%
  dplyr::mutate(Diff_First=First_18Pop_Pct-Prev_First,
                Diff_Complete=Complete_18Pop_Pct-Prev_Complete)
### Subset to CDC Midwest region, and start 4 weeks prior to first lottery (Week 15): 
Vax.dat <- Vax_wk_2021 %>% 
  dplyr::filter(Cluster %in% c("IA","IL","IN","KS","MI","MN","MO","ND",
                             "NE","OH","SD","WI"),
                Period >= 15)

## Get unique periods, clusters, and start times:
Periods <- unique(Vax.dat$Period)
OrderedPds <- Periods[order(Periods)]
J <- length(OrderedPds)

## Prep Outcome Data in appropriate order and get StartTimes and N:
Ord_Data <- Vax.dat %>% dplyr::rename(StartPd=lott_week) %>%
  arrange(StartPd,Cluster,Period)
StartTimes <- Ord_Data %>% dplyr::select(Cluster,StartPd) %>%
  distinct()
N <- length(StartTimes$Cluster)
Obs_Y <- matrix(data=c(Ord_Data$First_18Pop_Pct, 
                       Ord_Data$Complete_18Pop_Pct,
                       Ord_Data$Diff_First,
                       Ord_Data$Diff_Complete), ncol=4)
colnames(Obs_Y) <- c("First_18Pop_Pct","Complete_18Pop_Pct",
                     "First_Pct_Diff","Complete_Pct_Diff")

## Generate A matrix:
Amat <- gen_A(N,J)

## Generate Theta object to see Theta values:
Theta2 <- gen_Theta(gen_js(StartTimes$Cluster,
                        StartTimes$StartPd,
                        OrderedPds),
                 Assumption=2)

## Run Solver for different assumption settings:
SO2 <- Solve_Assumption(Amat, StartTimes, OrderedPds,
                        Assumption=2,
                        v.Mat=cbind(Avg=rep(1/26,26),
                                    D.1=create_V(26, c(1,6,10,19)),
                                    D.2=create_V(26, c(2,8,13,23)),
                                    D.3=create_V(26, c(3, 11, 16)),
                                    D.4=create_V(26, c(4, 14, 20)),
                                    D.12=create_V(26, c(1, 2, 6, 8, 10, 13, 19, 23)),
                                    D.1234=create_V(26, c(1:4,6,8,11,14,10,13,16,20)),
                                    D.234=create_V(26, c(2:4, 8, 11, 14, 13, 16, 20)),
                                    OH=create_V(26, c(1:5, 7, 9, 12, 15, 18, 22, 26)),
                                    IL=create_V(26, c(6, 8, 11, 14, 17, 21, 25)),
                                    MI=create_V(26, c(10, 13, 16, 20, 24)),
                                    MO=create_V(26, c(19, 23)),
                                    OH2=create_V(26, 2),
                                    IL2=create_V(26, 8),
                                    MI2=create_V(26, 13),
                                    MO2=create_V(26, 23),
                                    Group=1/4*(create_V(26, c(1:5, 7, 9, 12, 15, 18, 22, 26))+
                                                 create_V(26, c(6, 8, 11, 14, 17, 21, 25))+
                                                 create_V(26, c(10, 13, 16, 20, 24))+
                                                 create_V(26, c(19, 23)))),
                        save_loc="../int_large/",
                        save_prefix="vax-solve-a_")
SO3 <- Solve_Assumption(Amat, StartTimes, OrderedPds,
                        Assumption=3,
                        v.Mat=cbind(D.1=create_V(12, 1),
                                    D.2=create_V(12, 2),
                                    D.12=create_V(12, 1:2),
                                    D.1234=create_V(12, 1:4),
                                    D.234=create_V(12, 2:4),
                                    Avg=create_V(12, 1:12)),
                        save_loc="../int_large/",
                        save_prefix="vax-solve-a_")
SO4 <- Solve_Assumption(Amat, StartTimes, OrderedPds,
                        Assumption=4,
                        v.Mat=cbind(Avg=create_V(12, 1:12),
                                    T.30=create_V(12, 12)),
                        save_loc="../int_large/",
                        save_prefix="vax-solve-a_")
SO5 <- Solve_Assumption(Amat, StartTimes, OrderedPds,
                        Assumption=5,
                        v.Mat=1,
                        save_loc="../int_large/",
                        save_prefix="vax-solve-a_")

## Run variance minimizer for different settings:
set.seed(1011)
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
                             save_prefix="vax-mv-a_"))
}

#### Estimate AR correlation on unused non-lottery states
ar.dat <- Vax_wk_2021 %>% 
  dplyr::filter(!(Cluster %in% c("IA","IL","IN","KS","MI","MN","MO","ND","NE","OH","SD","WI")),
                lottery==0,
                Period >= 15)
require(geepack)
GEE <- geepack::geeglm(First_18Pop_Pct~as.factor(Period)+as.factor(Cluster),
                       data=ar.dat,
                       family=gaussian,
                       id=Cluster,
                       corstr="ar1")

### AR(1) (rho = 0.95 estimated from above; First_Diff has 0.34 instead):
for (i in 2:5) {
  assign(x=paste0("MVOut_",i,"_AR1_0_95"),
         value=MV_Assumption(SolveOut=get(paste0("SO",i)),
                             Assumption=i,
                             Sigma=create_Sigma_AR1(rho=0.95,N=N,J=J),
                             SigmaName="AR1_0_95",
                             Observations=Obs_Y,
                             Permutations=1000,
                             save_loc="int/",
                             save_prefix="vax-mv-a_"))
}

## Import Results:
Assns <- 2:5
SigmaNames <- c("Ind","AR1_0_95")

for (j in SigmaNames) {
  for (i in Assns) {
    load(file=paste0("int/vax-mv-a_",i,"_",j,".Rda"))
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
  FirstEsts <- NULL
  CompleteEsts <- NULL
  Diff1Ests <- NULL
  DiffCEsts <- NULL
  for (j in SigmaNames) {
    FirstEsts <- cbind(FirstEsts,(get(paste0("MVOut_",i,"_",j))[["Estimates"]])[,"First_18Pop_Pct"])
    CompleteEsts <- cbind(CompleteEsts,exp((get(paste0("MVOut_",i,"_",j))[["Estimates"]])[,"Complete_18Pop_Pct"]))
    Diff1Ests <- cbind(Diff1Ests,(get(paste0("MVOut_",i,"_",j))[["Estimates"]])[,"First_Pct_Diff"])
    DiffCEsts <- cbind(DiffCEsts,exp((get(paste0("MVOut_",i,"_",j))[["Estimates"]])[,"Complete_Pct_Diff"]))
  }
  colnames(FirstEsts) <- SigmaNames
  colnames(CompleteEsts) <- SigmaNames
  colnames(Diff1Ests) <- SigmaNames
  colnames(DiffCEsts) <- SigmaNames
  print("First Dose Percentage Estimates:")
  print(FirstEsts)
  print("Complete Series Percentage Estimates:")
  print(CompleteEsts)
  print("First Dose Percentage Difference Estimates:")
  print(Diff1Ests)
  print("Complete Series Percentage Difference Estimates:")
  print(DiffCEsts)
}

for (i in Assns) {
  print(paste0("P-Values for Assumption ",i))
  for (j in SigmaNames) {
    print(j)
    print((get(paste0("MVOut_",i,"_",j)))[["P_Values"]])
  }
}

## Print Observation Weight Heatmaps:
### To create various heat maps, add rows with 
### different values of i (Assumption Setting),
### j (Variance setting), and Estimators (estimator)
Map_Settings <- tibble(i=c(rep(2,6),5),
                       j=rep("AR1_0_95",7),
                       Estimators=c("Avg","D.2","D.1234","D.234","OH","Group",1),
                       Est_labs=c("Overall ATT", "Second-Period Effect",
                                  "First Four Weeks Effect", "Weeks 2-4 Effect", "Ohio",
                                  "State-Averaged","Overall, Assumption S5"))
for (row in 1:(dim(Map_Settings)[1])) {
  Weights <- (get(paste0("MVOut_",Map_Settings[row,] %>% pull("i"),"_",
                         Map_Settings[row,] %>% pull("j")))[["MV"]])[["Obs.weights"]]
  if (is.null(colnames(Weights))) {
    colnames(Weights) <- as.character(1:(dim(Weights)[2]))
  }
  Obs.weight.dat <- tibble(x=rep(1:J, times=N), y=rep(1:N, each=J),
                           Value=Weights[,Map_Settings[row,] %>% pull("Estimators")])
  ggsave(filename=paste0("figs/Vax-Weights_Heatmap_",Map_Settings[row,"i"],"_",
                         Map_Settings[row,] %>% pull("j"),"_",
                         Map_Settings[row,] %>% pull("Estimators"),".png"),
         plot=ggplot(data=Obs.weight.dat, mapping=aes(x=x, y=y, fill=Value)) +
           geom_tile() + theme_bw() +
           coord_cartesian(xlim=c(0.5,J+0.5), ylim=c(N+0.5,0.5), clip="off", expand=FALSE) +
           scale_y_reverse(breaks=1:N, 
                           labels=StartTimes$Cluster,
                           minor_breaks=NULL) +
           scale_x_continuous(breaks=1:J, 
                              labels=OrderedPds,
                              minor_breaks=NULL) +
           scale_fill_gradient2(low="#542788",high="#b35806") +
           labs(x="MMWR Week (2021)", y="State", fill="Weight",
                title=paste0("Observation Weights, ", 
                             Map_Settings[row,] %>% pull("Est_labs"))),
         width=6, height=4, units="in", dpi=600)
}
  

## Comparisons to other methods:
### Get comparison estimates:
DFT <- SO5$DFT
Comp_wts <- Comp_Ests_Weights(DFT_obj=DFT, Amat=Amat,
                              estimator=c("TW","CS","SA","CH","CO","NP"))
Comp_ests <- t(as.matrix(Comp_wts$Obs.weights)) %*% Obs_Y
Comp_ests

### Get comparison perm. p-values
set.seed(12959)
Comp_perms <- replicate(n=1000,
                        expr=Permute_obs(Observations=Obs_Y, 
                                         N=DFT$N, J=DFT$J, 
                                         Obs.weights=Comp_wts$Obs.weights))
Comp_perms2 <- simplify2array(apply(Comp_perms, 3, 
                                    FUN=function(x) abs(x) >= abs(Comp_ests), 
                                    simplify=FALSE))
Comp_pvals <- apply(Comp_perms2, c(1,2), mean)
Comparisons=list(Estimates=Comp_ests, P_Values=Comp_pvals)

### Save comparisons:
save(Comparisons, file="int/Vax-Comp-Ests.Rda")

## Check against existing packages for staggered adoption methods:
### Packages:
require(did) ## For CS
require(fixest) ## For SA
require(DIDmultiplegt) ## For CH

#### Data Prep:
Cl.labs <- StartTimes %>%
  mutate(Cl.Num=1:N)
Vax.dat.2 <- Ord_Data %>%
  left_join(Cl.labs %>% dplyr::select(-c(StartPd)), by=join_by(Cluster)) %>%
  dplyr::mutate(Interv=Interv==1)

### TWFE:
TW <- lm(First_18Pop_Pct~Interv+factor(Period)+Cluster, data=Vax.dat.2)
coef(TW)["IntervTRUE"]

### Callaway and Sant'Anna (2021):
CS_gt <- att_gt(yname="First_18Pop_Pct",
                tname="Period",
                idname="Cl.Num",
                gname="StartPd",
                data=Vax.dat.2,
                panel=TRUE,
                control_group="notyettreated")
ggdid(CS_gt)
CS_gt
aggte(CS_gt, type="simple")
aggte(CS_gt, type="dynamic")
aggte(CS_gt, type="group")
aggte(CS_gt, type="calendar")


### Sun and Abraham (2021):
SA <- feols(First_18Pop_Pct~sunab(cohort=StartPd,
                                  period=Period, 
                                  att=TRUE) | Cluster + Period,
            data=Vax.dat.2)
summary(SA)

### de Chaisemartin and d'Haultfoeuille (2020):
CH <- did_multiplegt(df=Vax.dat.2,
                     Y="First_18Pop_Pct",
                     G="Cluster",
                     T="Period",
                     D="Interv")
CH

