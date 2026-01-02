#######################################
###### File: Xpert-analysis.R #########
###### Lee Kennedy-Shaffer ############
#######################################

require(tidyverse)
require(lme4)
set.seed(413354)

source("R/Full_Analysis.R")
source("R/CompEsts.R")

## Read in (simulated) data from data folder:
load("data/Xpert-data-sim.Rda")

## Get unique periods, start times, J:
Periods <- unique(xpert.dat$Period)
OrderedPds <- Periods[order(Periods)]
J <- length(OrderedPds)

## Create StartTimes DF and get N:
StartTimes <- xpert.dat %>% dplyr::filter(Interv==1) %>%
  group_by(Cluster) %>% dplyr::summarize(StartPd=min(Period)) %>%
  dplyr::arrange(StartPd,Cluster)
N <- nrow(StartTimes)

## Prep Outcome Data in appropriate order:
Ord_Data <- xpert.dat %>% left_join(StartTimes, by="Cluster") %>%
  dplyr::arrange(StartPd,Cluster,Period)
Obs_Y <- matrix(data=c(Ord_Data$Outcome, Ord_Data$logOdds), ncol=2)
colnames(Obs_Y) <- c("Probability","Log Odds")

## Creating Schematic Figure:
ggsave(filename="figs/Xpert_Schematic.eps",
       plot=ggplot(data=Ord_Data %>%
                     dplyr::left_join(StartTimes %>% dplyr::mutate(ClusterNum=1:N),
                                                      by=join_by(Cluster,StartPd)) %>%
                     dplyr::mutate(Treatment=factor(Interv,
                                                    levels=c(0,1),
                                                    labels=c("Control","Intervention"))),
                   mapping=aes(x=Period, y=ClusterNum, fill=Treatment)) +
         geom_tile(color="grey80", lty=1) + theme_bw() +
         coord_cartesian(xlim=c(0.5,J+0.5), ylim=c(N+0.5,0.5),
                         clip="off", expand=FALSE) +
         scale_y_reverse(breaks=1:N, minor_breaks=NULL) +
         scale_x_continuous(breaks=1:J, minor_breaks=NULL) +
         scale_fill_manual(values=c("white","grey20")) +
         labs(x="Period (Month of Study)", y="Cluster", fill="Treatment",
              title="Trial Schematic for Example SWT"),
       width=6, height=4, units="in")

## Run Solver for different assumption settings:
SO2 <- Solve_Assumption(StartTimes,
                        OrderedPds,
                        Assumption=2,
                        v.Mat=cbind(Avg=c(rep(1/28,28)),
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
                                                rep(0,7))))

SO3 <- Solve_Assumption(StartTimes,
                        OrderedPds,
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

SO4 <- Solve_Assumption(StartTimes,
                        OrderedPds,
                        Assumption=4,
                        v.Mat=cbind(AvgEx8=c(rep(1/6,6),0),
                             T.2=c(1,rep(0,6)),
                             T.3=c(0,1,rep(0,5)),
                             T.4=c(0,0,1,rep(0,4)),
                             T.5=c(0,0,0,1,rep(0,3)),
                             T.6=c(rep(0,4),1,0,0),
                             T.7=c(rep(0,5),1,0),
                             T.8=c(rep(0,6),1),
                             Middle=c(0,1/3,1/3,1/3,0,0,0)))

SO5 <- Solve_Assumption(StartTimes,
                        OrderedPds,
                        Assumption=5,
                        v.Mat=1)

## Run variance minimizer for different settings:
Assns <- 2:5
SigmaNames <- c("Ind","CS_0_003","CS_0_333","AR1_0_012")

for (SigName in SigmaNames) {
  if (SigName=="Ind") {
    Sig <- create_Sigma_Ind(N=N, J=J)
  } else if (SigName=="CS_0_003") {
    ### Using rho = 0.003 from Thompson et al. 2018
    Sig <- create_Sigma_CS(rho=0.003, N=N, J=J)
  } else if (SigName=="CS_0_333") {
    ### Alternative version to test robustness to input ICC
    Sig <- create_Sigma_CS(rho=1/3, N=N, J=J)
  } else if (SigName=="AR1_0_012") {
    ### (rho = 0.012 gives average ICC within a cluster ~0.003):
    Sig <- create_Sigma_AR1(rho=0.012, N=N, J=J)
  }
  for (i in Assns) {
    assign(x=paste("MVOut",i,SigName, sep="_"),
           value=MV_Assumption(SolveOut=get(paste0("SO",i)),
                               Assumption=i,
                               Sigma=Sig,
                               SigmaName=SigName,
                               Observations=Obs_Y,
                               Permutations=1000))
  }
}

### Independence:
# for (i in 2:5) {
#   assign(x=paste0("MVOut_",i,"_Ind"),
#          value=MV_Assumption(SolveOut=get(paste0("SO",i)),
#                              Assumption=i,
#                              Sigma=create_Sigma_Ind(N=N,J=J),
#                              SigmaName="Ind",
#                              Observations=Obs_Y,
#                              Permutations=1000 #,
#                              # save_loc="int_new/",
#                              # save_prefix="xpert-mv-a_"
#                              ))
# }

## Export/Import Results:

for (j in SigmaNames) {
  for (i in Assns) {
    save(list=paste("MVOut",i,j, sep="_"),
         file=paste0("int/xpert-mv-a_",i,"_",j,".Rda"))
  }
}

# for (j in SigmaNames) {
#   for (i in Assns) {
#     load(file=paste0("int/xpert-mv-a_",i,"_",j,".Rda"))
#   }
# }

## Summarize Results:
for (j in SigmaNames) {
  print(paste0("Variance: ",j))
  for (i in Assns) {
    print(paste0("Assumption: ",i))
    print((get(paste0("MVOut_",i,"_",j))[["MV"]])[["Variance"]])
  }
}

for (i in Assns) {
  print(paste0("Assumption: ",i))
  ProbEsts <- NULL
  OREsts <- NULL
  for (j in SigmaNames) {
    ProbEsts <- cbind(ProbEsts,(get(paste0("MVOut_",i,"_",j))[["Estimates"]])[,"Probability"])
    OREsts <- cbind(OREsts,exp((get(paste0("MVOut_",i,"_",j))[["Estimates"]])[,"Log Odds"]))
  }
  colnames(ProbEsts) <- SigmaNames
  colnames(OREsts) <- SigmaNames
  print("Risk Difference Estimates:")
  print(ProbEsts)
  print("Odds Ratio:")
  print(OREsts)
}


### Inference:
for (i in Assns) {
  print(paste0("P-Values for Assumption: ",i))
  for (j in SigmaNames) {
    print(j)
    print((get(paste0("MVOut_",i,"_",j)))[["P_Values"]])
  }
}

### Observation Weight Heatmaps:
### To create various heat maps, add rows with
### different values of i (Assumption Setting),
### j (Variance setting), and Estimators (estimator)
Map_Settings <- tibble(i=c(5,4,3,2,3,2,rep(4,6),rep(2,6)),
                       j=rep("CS_0_003",18),
                       Estimators=c("1","AvgEx8","Avg","AvgExT8",
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
    Weights <- (get(paste0("MVOut_",Map_Settings[row,] %>% pull("i"),"_",
                           Map_Settings[row,] %>% pull("j")))[["MV"]])[["Obs.weights"]]
    if (is.null(colnames(Weights))) {
      colnames(Weights) <- as.character(1:(dim(Weights)[2]))
    }
    Obs.weight.dat <- tibble(x=rep(1:J, times=N), y=rep(1:N, each=J),
                             Value=Weights[,Map_Settings[row,] %>% pull("Estimators")])
    ggsave(filename=paste0("figs/Xpert-Weights_Heatmap_",Map_Settings[row,"i"],"_",
                           Map_Settings[row,] %>% pull("j"),"_",
                           Map_Settings[row,] %>% pull("Estimators"),".eps"),
           plot=ggplot(data=Obs.weight.dat, mapping=aes(x=x, y=y, fill=Value)) +
             geom_tile(color="grey80", lty=1) + theme_bw() +
             geom_text(aes(label=format(round(Value, digits=3), nsmall=3))) +
             coord_cartesian(xlim=c(0.5,J+0.5), ylim=c(N+0.5,0.5), clip="off", expand=FALSE) +
             scale_y_reverse(breaks=1:N, minor_breaks=NULL) +
             scale_x_continuous(breaks=1:J, minor_breaks=NULL) +
             scale_fill_gradient2(low="#7938C6",high="#DE6D07") +
             labs(x="Period (Month of Study)", y="Cluster", fill="Weight",
                  title=paste0("Observation Weights: ",
                               Map_Settings[row,] %>% pull("Est_labs"))),
           width=6, height=4, units="in")
}

## Comparisons to other methods:
### Get comparison estimates for methods with known weights:
Comp_wts <- Comp_Ests_Weights(ADFT_obj=SO5$ADFT,
                              estimator=c("TW","CS","SA","CH","CO","NP"))
Comp_ests <- t(as.matrix(Comp_wts$Obs.weights)) %*% Obs_Y
Comp_ests

### Get comparison estimates and perm. p-values for CLWP/CLWPA:
xpert.dat.long <- NULL
for (i in 1:(dim(xpert.dat)[1])) {
  add <- xpert.dat[i,,drop=FALSE] %>% dplyr::select(Interv,Period,Cluster) %>%
    cross_join(tibble(Outcome=c(rep(1,xpert.dat[i,"Events"]),
                                rep(0,xpert.dat[i,"Indivs"]-xpert.dat[i,"Events"]))))
  xpert.dat.long <- xpert.dat.long %>% bind_rows(add)
}

GetCLWPs <- function(data,start_vals) {
  Fit_Lin <- CLWP_fit(my.data=data %>% rename(Y.ij.bar=Outcome),
                      start_theta=start_vals["theta_lin"],
                      start_sigma=start_vals["sigma_lin"],
                      N=N)
  Fit_A_Lin <- CLWPA_fit(my.data=data %>% rename(Y.ij.bar=Outcome),
                         start_theta=start_vals["theta_lin"],
                         start_sigma=start_vals["sigma_lin"],
                         N=N)
  Fit_Log <- CLWP_fit(my.data=data %>% rename(Y.ij.bar=logOdds),
                      start_theta=start_vals["theta_log"],
                      start_sigma=start_vals["sigma_log"],
                      N=N)
  Fit_A_Log <- CLWPA_fit(my.data=data %>% rename(Y.ij.bar=logOdds),
                         start_theta=start_vals["theta_log"],
                         start_sigma=start_vals["sigma_log"],
                         N=N)
  Ests <- matrix(data=unname(c(Fit_Lin["CLWP_est.theta"],
                               Fit_A_Lin["CLWPA_est.theta"],
                               Fit_Log["CLWP_est.theta"],
                               Fit_A_Log["CLWPA_est.theta"])),
                 nrow=2,ncol=2,byrow=FALSE)
  rownames(Ests) <- c("CLWP","CLWPA")
  colnames(Ests) <- c("Probability","Log Odds")
  return(Ests)
}

MEM.start <- lmer(Outcome~Interv+factor(Period)+(1|Cluster),
                  data=xpert.dat.long)
start_theta_lin <- unname(fixef(MEM.start)["Interv"])
start_sigma_lin <- sqrt(as.data.frame(VarCorr(MEM.start))[2,"vcov"]*2)
MEM.start.log <- glmer(Outcome~Interv+factor(Period)+(1|Cluster),
                       data=xpert.dat.long,
                       family=binomial)
start_theta_log <- unname(fixef(MEM.start.log)["Interv"])
start_sigma_log <- sd(summary(MEM.start.log)$residuals)*sqrt(2)
start_vals <- c(start_theta_lin,start_sigma_lin,
                start_theta_log,start_sigma_log)
names(start_vals) <- c("theta_lin","sigma_lin",
                       "theta_log","sigma_log")

Comp_ests <- rbind(Comp_ests,
                   GetCLWPs(data=xpert.dat, start_vals=start_vals))

### Get comparison perm. p-values:
set.seed(7446)
SinglePerm <- function() {
  Perm_out <- Permute_obs(Observations=Obs_Y,
                          N=SO5$ADFT$N, J=SO5$ADFT$J,
                          Obs.weights=Comp_wts$Obs.weights)
  return(rbind(Perm_out$Ests,
               GetCLWPs(data=xpert.dat %>%
                          dplyr::select(Interv,Period,Cluster) %>%
                          bind_cols(Perm_out$Obs) %>%
                          rename(Outcome=Probability,
                                 logOdds=`Log Odds`),
                        start_vals=start_vals)))
}
Comp_perms <- replicate(n=1000, expr=SinglePerm())
Comp_perms2 <- simplify2array(apply(Comp_perms, 3,
                                    FUN=function(x) abs(x) >= abs(Comp_ests), simplify=FALSE))
Comp_pvals <- apply(Comp_perms2, c(1,2), mean)
Comparisons=list(Estimates=Comp_ests, P_Values=Comp_pvals)

### Save comparisons:
save(Comparisons, file="int/xpert-Comp-Ests.Rda")


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

### TWFE:
TW <- lm(Outcome~Interv+factor(Period)+ClusterF, data=xpert.dat.2)
coef(TW)["Interv"]

### Callaway and Sant'Anna (2021):
CS_gt <- att_gt(yname="Outcome",
                tname="Period",
                idname="Cluster",
                gname="StartPd",
                data=xpert.dat.2,
                panel=TRUE,
                control_group="notyettreated")
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
CH <- did_multiplegt(mode="old",
                     df=xpert.dat,
                     Y="Outcome",
                     G="Cluster",
                     T="Period",
                     D="Interv")
CH

## MS Results:

Tbl2 <- tibble(Assumption=5:2,
               Estimator=c("1","AvgEx8","Avg","AvgExT8"),
               OR_Est=exp(c(MVOut_5_CS_0_003$Estimates[1,"Log Odds"],
                            MVOut_4_CS_0_003$Estimates["AvgEx8","Log Odds"],
                            MVOut_3_CS_0_003$Estimates["Avg","Log Odds"],
                            MVOut_2_CS_0_003$Estimates["AvgExT8","Log Odds"])),
               OR_P=c(MVOut_5_CS_0_003$P_Values[1,"Log Odds"],
                      MVOut_4_CS_0_003$P_Values["AvgEx8","Log Odds"],
                      MVOut_3_CS_0_003$P_Values["Avg","Log Odds"],
                      MVOut_2_CS_0_003$P_Values["AvgExT8","Log Odds"]),
               RD_Est=100*c(MVOut_5_CS_0_003$Estimates[1,"Probability"],
                            MVOut_4_CS_0_003$Estimates["AvgEx8","Probability"],
                            MVOut_3_CS_0_003$Estimates["Avg","Probability"],
                            MVOut_2_CS_0_003$Estimates["AvgExT8","Probability"]),
               RD_P=c(MVOut_5_CS_0_003$P_Values[1,"Probability"],
                      MVOut_4_CS_0_003$P_Values["AvgEx8","Probability"],
                      MVOut_3_CS_0_003$P_Values["Avg","Probability"],
                      MVOut_2_CS_0_003$P_Values["AvgExT8","Probability"])
)

Vars <- c(MVOut_5_CS_0_003$MV$Variance,
          MVOut_4_CS_0_003$MV$Variance[1,"AvgEx8"],
          MVOut_3_CS_0_003$MV$Variance[1,"Avg"],
          MVOut_2_CS_0_003$MV$Variance[1,"AvgExT8"])
Rel_Effs <- c(Vars[2]/Vars[1], Vars[3]/Vars[1], Vars[4]/Vars[1])

Tbl3 <- tibble(Month=2:7,
               OR_S4_Est=exp(c(MVOut_4_CS_0_003$Estimates["T.2","Log Odds"],
                               MVOut_4_CS_0_003$Estimates["T.3","Log Odds"],
                               MVOut_4_CS_0_003$Estimates["T.4","Log Odds"],
                               MVOut_4_CS_0_003$Estimates["T.5","Log Odds"],
                               MVOut_4_CS_0_003$Estimates["T.6","Log Odds"],
                               MVOut_4_CS_0_003$Estimates["T.7","Log Odds"])),
               OR_S4_P=c(MVOut_4_CS_0_003$P_Values["T.2","Log Odds"],
                         MVOut_4_CS_0_003$P_Values["T.3","Log Odds"],
                         MVOut_4_CS_0_003$P_Values["T.4","Log Odds"],
                         MVOut_4_CS_0_003$P_Values["T.5","Log Odds"],
                         MVOut_4_CS_0_003$P_Values["T.6","Log Odds"],
                         MVOut_4_CS_0_003$P_Values["T.7","Log Odds"]),
               OR_S2_Est=exp(c(MVOut_2_CS_0_003$Estimates["T.2","Log Odds"],
                               MVOut_2_CS_0_003$Estimates["T.3","Log Odds"],
                               MVOut_2_CS_0_003$Estimates["T.4","Log Odds"],
                               MVOut_2_CS_0_003$Estimates["T.5","Log Odds"],
                               MVOut_2_CS_0_003$Estimates["T.6","Log Odds"],
                               MVOut_2_CS_0_003$Estimates["T.7","Log Odds"])),
               OR_S2_P=c(MVOut_2_CS_0_003$P_Values["T.2","Log Odds"],
                         MVOut_2_CS_0_003$P_Values["T.3","Log Odds"],
                         MVOut_2_CS_0_003$P_Values["T.4","Log Odds"],
                         MVOut_2_CS_0_003$P_Values["T.5","Log Odds"],
                         MVOut_2_CS_0_003$P_Values["T.6","Log Odds"],
                         MVOut_2_CS_0_003$P_Values["T.7","Log Odds"]),
               RD_S4_Est=100*c(MVOut_4_CS_0_003$Estimates["T.2","Probability"],
                               MVOut_4_CS_0_003$Estimates["T.3","Probability"],
                               MVOut_4_CS_0_003$Estimates["T.4","Probability"],
                               MVOut_4_CS_0_003$Estimates["T.5","Probability"],
                               MVOut_4_CS_0_003$Estimates["T.6","Probability"],
                               MVOut_4_CS_0_003$Estimates["T.7","Probability"]),
               RD_S4_P=c(MVOut_4_CS_0_003$P_Values["T.2","Probability"],
                         MVOut_4_CS_0_003$P_Values["T.3","Probability"],
                         MVOut_4_CS_0_003$P_Values["T.4","Probability"],
                         MVOut_4_CS_0_003$P_Values["T.5","Probability"],
                         MVOut_4_CS_0_003$P_Values["T.6","Probability"],
                         MVOut_4_CS_0_003$P_Values["T.7","Probability"]),
               RD_S2_Est=100*c(MVOut_2_CS_0_003$Estimates["T.2","Probability"],
                               MVOut_2_CS_0_003$Estimates["T.3","Probability"],
                               MVOut_2_CS_0_003$Estimates["T.4","Probability"],
                               MVOut_2_CS_0_003$Estimates["T.5","Probability"],
                               MVOut_2_CS_0_003$Estimates["T.6","Probability"],
                               MVOut_2_CS_0_003$Estimates["T.7","Probability"]),
               RD_S2_P=c(MVOut_2_CS_0_003$P_Values["T.2","Probability"],
                         MVOut_2_CS_0_003$P_Values["T.3","Probability"],
                         MVOut_2_CS_0_003$P_Values["T.4","Probability"],
                         MVOut_2_CS_0_003$P_Values["T.5","Probability"],
                         MVOut_2_CS_0_003$P_Values["T.6","Probability"],
                         MVOut_2_CS_0_003$P_Values["T.7","Probability"])
)

TblWA3 <- tibble(Assumption=5:2,
                 Estimator=c("1","AvgEx8","Avg","AvgExT8"),
                 Ind_Est=exp(c(MVOut_5_Ind$Estimates[1,"Log Odds"],
                               MVOut_4_Ind$Estimates["AvgEx8","Log Odds"],
                               MVOut_3_Ind$Estimates["Avg","Log Odds"],
                               MVOut_2_Ind$Estimates["AvgExT8","Log Odds"])),
                 Ind_P=c(MVOut_5_Ind$P_Values[1,"Log Odds"],
                         MVOut_4_Ind$P_Values["AvgEx8","Log Odds"],
                         MVOut_3_Ind$P_Values["Avg","Log Odds"],
                         MVOut_2_Ind$P_Values["AvgExT8","Log Odds"]),
                 CS_Est=exp(c(MVOut_5_CS_0_003$Estimates[1,"Log Odds"],
                              MVOut_4_CS_0_003$Estimates["AvgEx8","Log Odds"],
                              MVOut_3_CS_0_003$Estimates["Avg","Log Odds"],
                              MVOut_2_CS_0_003$Estimates["AvgExT8","Log Odds"])),
                 CS_P=c(MVOut_5_CS_0_003$P_Values[1,"Log Odds"],
                        MVOut_4_CS_0_003$P_Values["AvgEx8","Log Odds"],
                        MVOut_3_CS_0_003$P_Values["Avg","Log Odds"],
                        MVOut_2_CS_0_003$P_Values["AvgExT8","Log Odds"]),
                 AR1_Est=exp(c(MVOut_5_AR1_0_012$Estimates[1,"Log Odds"],
                               MVOut_4_AR1_0_012$Estimates["AvgEx8","Log Odds"],
                               MVOut_3_AR1_0_012$Estimates["Avg","Log Odds"],
                               MVOut_2_AR1_0_012$Estimates["AvgExT8","Log Odds"])),
                 AR1_P=c(MVOut_5_AR1_0_012$P_Values[1,"Log Odds"],
                         MVOut_4_AR1_0_012$P_Values["AvgEx8","Log Odds"],
                         MVOut_3_AR1_0_012$P_Values["Avg","Log Odds"],
                         MVOut_2_AR1_0_012$P_Values["AvgExT8","Log Odds"])
)

TblWA4 <- tibble(Assumption=5:2,
                 Estimator=c("1","AvgEx8","Avg","AvgExT8"),
                 Ind_Est=100*c(MVOut_5_Ind$Estimates[1,"Probability"],
                               MVOut_4_Ind$Estimates["AvgEx8","Probability"],
                               MVOut_3_Ind$Estimates["Avg","Probability"],
                               MVOut_2_Ind$Estimates["AvgExT8","Probability"]),
                 Ind_P=c(MVOut_5_Ind$P_Values[1,"Probability"],
                         MVOut_4_Ind$P_Values["AvgEx8","Probability"],
                         MVOut_3_Ind$P_Values["Avg","Probability"],
                         MVOut_2_Ind$P_Values["AvgExT8","Probability"]),
                 CS_Est=100*c(MVOut_5_CS_0_003$Estimates[1,"Probability"],
                              MVOut_4_CS_0_003$Estimates["AvgEx8","Probability"],
                              MVOut_3_CS_0_003$Estimates["Avg","Probability"],
                              MVOut_2_CS_0_003$Estimates["AvgExT8","Probability"]),
                 CS_P=c(MVOut_5_CS_0_003$P_Values[1,"Probability"],
                        MVOut_4_CS_0_003$P_Values["AvgEx8","Probability"],
                        MVOut_3_CS_0_003$P_Values["Avg","Probability"],
                        MVOut_2_CS_0_003$P_Values["AvgExT8","Probability"]),
                 AR1_Est=100*c(MVOut_5_AR1_0_012$Estimates[1,"Probability"],
                               MVOut_4_AR1_0_012$Estimates["AvgEx8","Probability"],
                               MVOut_3_AR1_0_012$Estimates["Avg","Probability"],
                               MVOut_2_AR1_0_012$Estimates["AvgExT8","Probability"]),
                 AR1_P=c(MVOut_5_AR1_0_012$P_Values[1,"Probability"],
                         MVOut_4_AR1_0_012$P_Values["AvgEx8","Probability"],
                         MVOut_3_AR1_0_012$P_Values["Avg","Probability"],
                         MVOut_2_AR1_0_012$P_Values["AvgExT8","Probability"])
)

TblWA5 <- tibble(Method=rownames(Comparisons$Estimates),
                 Estimate=100*Comparisons$Estimates[,"Probability"],
                 P_Value=Comparisons$P_Values[,"Probability"],
                 Estimator_GD=c("S5","S2_AvgExT8","S3_Avg",
                                "S2_Group","S4_AvgEx8",
                                "S2_AvgExT8","S2_D.1","S2_D.1",
                                "S3_D.1","S5",rep(NA_character_,5)),
                 Estimate_GD=100*c(MVOut_5_CS_0_003$Estimates[1,"Probability"],
                                   MVOut_2_CS_0_003$Estimates["AvgExT8","Probability"],
                                   MVOut_3_CS_0_003$Estimates["Avg","Probability"],
                                   MVOut_2_CS_0_003$Estimates["Group","Probability"],
                                   MVOut_4_CS_0_003$Estimates["AvgEx8","Probability"],
                                   MVOut_2_CS_0_003$Estimates["AvgExT8","Probability"],
                                   MVOut_2_CS_0_003$Estimates["D.1","Probability"],
                                   MVOut_2_CS_0_003$Estimates["D.1","Probability"],
                                   MVOut_3_CS_0_003$Estimates["D.1","Probability"],
                                   MVOut_5_CS_0_003$Estimates[1,"Probability"],
                                   rep(NA_real_,5)),
                 P_Value_GD=c(MVOut_5_CS_0_003$P_Values[1,"Probability"],
                              MVOut_2_CS_0_003$P_Values["AvgExT8","Probability"],
                              MVOut_3_CS_0_003$P_Values["Avg","Probability"],
                              MVOut_2_CS_0_003$P_Values["Group","Probability"],
                              MVOut_4_CS_0_003$P_Values["AvgEx8","Probability"],
                              MVOut_2_CS_0_003$P_Values["AvgExT8","Probability"],
                              MVOut_2_CS_0_003$P_Values["D.1","Probability"],
                              MVOut_2_CS_0_003$P_Values["D.1","Probability"],
                              MVOut_3_CS_0_003$P_Values["D.1","Probability"],
                              MVOut_5_CS_0_003$P_Values[1,"Probability"],
                              rep(NA_real_,5))
)

save(list=c("Tbl2","Vars","Rel_Effs","Tbl3",
            "TblWA3","TblWA4","TblWA5"),
     file="res/xpert_results.Rda")


