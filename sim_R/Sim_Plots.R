#######################################
###### File: Sim_Plots.R ##############
#######################################

library(tidyverse)
library(patchwork)

Colors <- c("#1b9e77","#d95f02","#7570b3","#e7298a")
Shapes <- c(15:18,8)

outdir <- "sim_res/figs/"

## Load simulation results:
load(file="sim_res/Simulation_Results.Rda")

SRP <- Sim_Results %>% dplyr::filter(Outcome=="Probability") %>%
  dplyr::mutate(Lower=`Mean Estimate`-`SD of Estimate`,
                Upper=`Mean Estimate`+`SD of Estimate`)

## Select Desired Estimators:
OverallSet <- tibble(Estimator=c("A5_Ind_Single","W_TW","CPI","W_CO.W_CO3",
                                 "A4_Ind_AvgEx8","CPI.T.TAvg","W_CS.W_calendar","A2_Ind_T.Avg","CPI.DT.TAvgExLast",
                                 "A3_Ind_Avg","CPI.D.DAvg","A3_Ind_AvgEx7","W_CS.W_dynamic","A2_Ind_D.Avg","CPI.DT.DAvg",
                                 "A2_Ind_Group","W_CS.W_group","CLWP","CLWPA",
                                 "A2_Ind_AvgExT8","W_CS.W_simple","W_SA.W_ATT","CPI.DT.DTAvgExLastT")) %>%
  mutate(`Estimator Number`=c(1:4,6:10,12:17,19:22,24:27),
         Type=c("GD","SA","ME","SA",
                "GD","ME","SA","GD","ME",
                "GD","ME","GD","SA","GD","ME",
                "GD","SA","CL","CL",
                "GD","SA","SA","ME"),
         Assumption=c("S5","S5","S5","S5",
                      "S4","S4","S4","S2","S2",
                      "S3","S3","S3","S3","S2","S2",
                      "S2","S2","S2","S2",
                      "S2","S2","S2","S2"),
         Name=c("GD_A5","TWFE","CPI_A5","CO3",
                "GD_A4","CPI_A4","CS_A4","GD_A2_T.Avg","CPI_A2_T.Avg",
                "GD_A3","CPI_A3","GD_A3_ExLast","CS_A3","GD_A2_D.Avg","CPI_A2_D.Avg",
                "GD_A2_group","CS_group","CLWP","CLWPA",
                "GD_A2_ATT","CS_ATT","SA_ATT","CPI_A2_AvgExLast"),
         Estimand=c(rep("Overall",4),
                    rep("Time Avg.",5),
                    rep("Exp. Avg.",6),
                    rep("Group Avg.",4),
                    rep("ATT",4)))
OverallSet_CS_0_333 <- tibble(Estimator=c("A5_CS_0_333_Single","W_TW","CPI","W_CO.W_CO3",
                                 "A4_CS_0_333_AvgEx8","CPI.T.TAvg","W_CS.W_calendar","A2_CS_0_333_T.Avg","CPI.DT.TAvgExLast",
                                 "A3_CS_0_333_Avg","CPI.D.DAvg","A3_CS_0_333_AvgEx7","W_CS.W_dynamic","A2_CS_0_333_D.Avg","CPI.DT.DAvg",
                                 "A2_CS_0_333_Group","W_CS.W_group","CLWP","CLWPA",
                                 "A2_CS_0_333_AvgExT8","W_CS.W_simple","W_SA.W_ATT","CPI.DT.DTAvgExLastT")) %>%
  mutate(`Estimator Number`=c(1:4,6:10,12:17,19:22,24:27),,
         Type=c("GD","SA","ME","SA",
                "GD","ME","SA","GD","ME",
                "GD","ME","GD","SA","GD","ME",
                "GD","SA","CL","CL",
                "GD","SA","SA","ME"),
         Assumption=c("S5","S5","S5","S5",
                      "S4","S4","S4","S2","S2",
                      "S3","S3","S3","S3","S2","S2",
                      "S2","S2","S2","S2",
                      "S2","S2","S2","S2"),
         Name=c("GD_A5","TWFE","CPI_A5","CO3",
                "GD_A4","CPI_A4","CS_A4","GD_A2_T.Avg","CPI_A2_T.Avg",
                "GD_A3","CPI_A3","GD_A3_ExLast","CS_A3","GD_A2_D.Avg","CPI_A2_D.Avg",
                "GD_A2_group","CS_group","CLWP","CLWPA",
                "GD_A2_ATT","CS_ATT","SA_ATT","CPI_A2_AvgExLast"),
         Estimand=c(rep("Overall",4),
                    rep("Time Avg.",5),
                    rep("Exp. Avg.",6),
                    rep("Group Avg.",4),
                    rep("ATT",4)))
TargetsSet <- tibble(Estimator=c("A4_Ind_T.3","CPI.T.Interv:PeriodF3","A2_Ind_T.3","CPI.DT.Pd3",
                                 "A3_Ind_D.2","CPI.D.Interv:DiffF2","A2_Ind_D.2","CPI.DT.Diff2",
                                 "A3_Ind_D.1","CPI.D.Interv:DiffF1","W_CO.W_CO2",
                                 "A2_Ind_D.1","CPI.DT.Diff1","W_CH.W_M","W_CO.W_CO1"
                                 )) %>%
  mutate(`Estimator Number`=c(1:4,6:9,11:17),
         Type=c("GD","ME","GD","ME",
                "GD","ME","GD","ME",
                "GD","ME","SA",
                "GD","ME","SA","SA"),
         Assumption=c("S4","S4","S2","S2",
                      "S3","S3","S2","S2",
                      "S3","S3","S2",
                      "S2","S2","S2","S2"),
         Name=c("GD_A4_j3","CPI_A4_j3","GD_A2_j3","CPI_A2_j3",
                "GD_A3_a2","CPI_A3_a2","GD_A2_a2","CPI_A2_a2",
                "GD_A3_a1","CPI_A3_a1","CO2",
                "GD_A2_a1","CPI_A2_a1","DCDH","CO1"),
         Estimand=c(rep("Time 3",4),
                    rep("Exp. Pd. 2",4),
                    rep("Exp. Pd. 1",7)))
TargetsSet_CS_0_333 <- tibble(Estimator=c("A4_CS_0_333_T.3","CPI.T.Interv:PeriodF3","A2_CS_0_333_T.3","CPI.DT.Pd3",
                                          "A3_CS_0_333_D.2","CPI.D.Interv:DiffF2","A2_CS_0_333_D.2","CPI.DT.Diff2",
                                          "A3_CS_0_333_D.1","CPI.D.Interv:DiffF1","W_CO.W_CO2",
                                          "A2_CS_0_333_D.1","CPI.DT.Diff1","W_CH.W_M","W_CO.W_CO1"
)) %>%
  mutate(`Estimator Number`=c(1:4,6:9,11:17),
         Type=c("GD","ME","GD","ME",
                "GD","ME","GD","ME",
                "GD","SA","ME",
                "GD","ME","SA","SA"),
         Assumption=c("S4","S4","S2","S2",
                      "S3","S3","S2","S2",
                      "S3","S3","S2",
                      "S2","S2","S2","S2"),
         Name=c("GD_A4_j3","CPI_A4_j3","GD_A2_j3","CPI_A2_j3",
                "GD_A3_a2","CPI_A3_a2","GD_A2_a2","CPI_A2_a2",
                "GD_A3_a1","CPI_A3_a1","CO2",
                "GD_A2_a1","CPI_A2_a1","DCDH","CO1"),
         Estimand=c(rep("Time 3",4),
                    rep("Exp. Pd. 2",4),
                    rep("Exp. Pd. 1",7)))

## Compile Desired Results:
Overall <- OverallSet %>% left_join(SRP %>% dplyr::select(-c(Outcome)),
                                     by=join_by(Estimator),
                                     multiple="all")
Overall_CS <- OverallSet_CS_0_333 %>% left_join(SRP %>% dplyr::select(-c(Outcome)),
                                                by=join_by(Estimator),
                                                multiple="all")
Target <- TargetsSet %>% left_join(SRP %>% dplyr::select(-c(Outcome)),
                                   by=join_by(Estimator),
                                   multiple="all")
Target_CS <- TargetsSet_CS_0_333 %>% left_join(SRP %>% dplyr::select(-c(Outcome)),
                                               by=join_by(Estimator),
                                               multiple="all")

Overall_Labels <- tibble(Name=c("Overall", "Calendar\nAveraged", "Exposure\nAveraged",
                                "Group\nAveraged", "Cluster-\nPeriod\nAveraged"),
                         Xval=c(2.5, 8, 14.5, 20.5, 25.5),
                         Yval=0.02)
Overall_Breaks <- c(0,5,11,18,23,28)

Target_Labels <- tibble(Name=c("Calendar\nPd. 3", "Exposure\nPd. 2", "Exposure\nPd. 1"),
                        Xval=c(2.5, 7.5, 14),
                        Yval=0.02)
Target_Breaks <- c(0,5,10,18)

### Functions to generate manuscript plots:

Power_Plots <- function(res_df, outname,
                        Labels, Breaks,
                        EstsR1, EstsR2, EstsR3) {
Power1 <-
  ggplot(res_df  %>% filter(Scenario==1,
                            Estimand %in% EstsR1),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=5, linetype="dashed", color="gray25") +
  geom_hline(yintercept=100*c(.05-3*sqrt(.05*.95/1000),
                              .05+3*sqrt(.05*.95/1000)),
             linetype="dotted", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Empirical Type I Error (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,10),
                     breaks=seq(0,10,by=2.5),
                     expand=c(0,0)) +

  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power2 <-
  ggplot(res_df  %>% filter(Scenario==2,
                            Estimand %in% EstsR1),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power3 <-
  ggplot(res_df  %>% filter(Scenario==3,
                            Estimand %in% EstsR1),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power4 <-
  ggplot(res_df  %>% filter(Scenario==4,
                            Estimand %in% EstsR2),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power5 <-
  ggplot(res_df  %>% filter(Scenario==5,
                            Estimand %in% EstsR2),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power6 <-
  ggplot(res_df  %>% filter(Scenario==6,
                            Estimand %in% EstsR2),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power7 <-
  ggplot(res_df  %>% filter(Scenario==7,
                             Estimand %in% EstsR3),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power8 <-
  ggplot(res_df  %>% filter(Scenario==8,
                             Estimand %in% EstsR3),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power9 <-
  ggplot(res_df  %>% filter(Scenario==9,
                             Estimand %in% EstsR3),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

ggsave(filename=paste0(outdir,paste0("Sim_Power_",outname,".eps")),
       plot=Power1 + guides(color="none", shape="none") + labs(x=NULL, title="A) Scenario 1") +
         Power2 + guides(color="none", shape="none") + labs(x=NULL, title="B) Scenario 2") +
         Power3 + guides(shape="none") +labs(x=NULL, title="C) Scenario 3") +
         Power4 + guides(color="none", shape="none") + labs(x=NULL, title="D) Scenario 4") +
         Power5 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="E) Scenario 5") +
         Power6 + guides(color="none") + labs(x=NULL, y=NULL, title="F) Scenario 6") +
         Power7 + guides(color="none", shape="none") + labs(title="G) Scenario 7") +
         Power8 + guides(color="none", shape="none") + labs(y=NULL, title="H) Scenario 8") +
         Power9 + guides(color="none", shape="none") + labs(y=NULL, title="I) Scenario 9") +
         plot_layout(nrow=3, ncol=3, byrow=TRUE, guides="collect") +
         theme(legend.position="bottom"),
       width=13, height=12, units="in")
}

Overall_Est_Plots <- function(res_df, outname,
                              Labels, Breaks) {
### Plots of Estimates and SDs:
Est1 <-
  ggplot(res_df %>% filter(Scenario==1),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=0, linetype="solid", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
  theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  scale_y_continuous(limits=c(-0.02,0.02),
                     breaks=seq(-0.02,0.02,by=0.005),
                     expand=c(0,0)) +
  scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
  geom_vline(xintercept=Breaks)

Est2 <-
  ggplot(res_df %>% filter(Scenario==2),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=-0.02, linetype="solid", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
  theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  scale_y_continuous(limits=c(-0.04,0.00),
                     breaks=seq(-0.04,0.00,by=0.005),
                     expand=c(0,0)) +
  scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
  geom_vline(xintercept=Breaks)

Est3 <-
  ggplot(res_df %>% filter(Scenario==3),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=-0.04, linetype="solid", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
  theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  scale_y_continuous(limits=c(-0.06, -0.02),
                     breaks=seq(-0.06, -0.02,by=0.005),
                     expand=c(0,0)) +
  scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
  geom_vline(xintercept=Breaks)

Est4 <-
  ggplot(res_df %>% filter(Scenario==4,
                           Estimand %in% c("Overall","Time Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_segment(y=-0.02, x=Breaks[2], xend=Breaks[3],
               linetype="solid", color="gray25") +
  geom_segment(y=0.005, x=Breaks[4], xend=Breaks[4]+2.5,
               linetype="solid", color="gray25") +
  geom_segment(y=-0.02, x=Breaks[5]-2.5, xend=Breaks[5],
               linetype="solid", color="gray25") +
  geom_segment(y=-0.0033333, x=Breaks[5], xend=Breaks[6],
               linetype="solid", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
  theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  scale_y_continuous(limits=c(-0.04,0.02),
                     breaks=seq(-0.04,0.02,by=0.01),
                     expand=c(0,0)) +
  scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
  geom_vline(xintercept=Breaks)

Est5 <-
  ggplot(res_df %>% filter(Scenario==5,
                           Estimand %in% c("Overall","Time Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_segment(y=-0.02, x=Breaks[2], xend=Breaks[3],
               linetype="solid", color="gray25") +
  geom_segment(y=0.005694444, x=Breaks[4], xend=Breaks[4]+2.5,
               linetype="solid", color="gray25") +
  geom_segment(y=-0.02, x=Breaks[5]-2.5, xend=Breaks[5],
               linetype="solid", color="gray25") +
  geom_segment(y=-0.001904762, x=Breaks[5], xend=Breaks[6],
               linetype="solid", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
  theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  scale_y_continuous(limits=c(-0.04,0.02),
                     breaks=seq(-0.04,0.02,by=0.01),
                     expand=c(0,0)) +
  scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
  geom_vline(xintercept=Breaks)

Est6 <-
  ggplot(res_df %>% filter(Scenario==6,
                           Estimand %in% c("Overall","Time Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_segment(y=-0.02, x=Breaks[2], xend=Breaks[3],
               linetype="solid", color="gray25") +
  geom_segment(y=-0.0105, x=Breaks[4], xend=Breaks[4]+2.5,
               linetype="solid", color="gray25") +
  geom_segment(y=-0.02, x=Breaks[5]-2.5, xend=Breaks[5],
               linetype="solid", color="gray25") +
  geom_segment(y=-0.01428571, x=Breaks[5], xend=Breaks[6],
               linetype="solid", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
  theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  scale_y_continuous(limits=c(-0.04,0.02),
                     breaks=seq(-0.04,0.02,by=0.01),
                     expand=c(0,0)) +
  scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
  geom_vline(xintercept=Breaks)

Est7 <-
  ggplot(res_df %>% filter(Scenario==7,
                           Estimand %in% c("Overall","Exp. Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_segment(y=-0.025, x=Breaks[3], xend=Breaks[3]+2.5,
               linetype="solid", color="gray25") +
  geom_segment(y=-0.0225, x=Breaks[3]+2.5, xend=Breaks[4],
               linetype="solid", color="gray25") +
  geom_segment(y=-.01625, x=Breaks[4], xend=Breaks[5],
               linetype="solid", color="gray25") +
  geom_segment(y=-0.01833333, x=Breaks[5], xend=Breaks[6],
               linetype="solid", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
  theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  scale_y_continuous(limits=c(-0.04,0.02),
                     breaks=seq(-0.04,0.02,by=0.01),
                     expand=c(0,0)) +
  scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
  geom_vline(xintercept=Breaks)

Est8 <-
  ggplot(res_df %>% filter(Scenario==8,
                           Estimand %in% c("Overall","Exp. Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_segment(y=-0.0214, x=Breaks[3], xend=Breaks[3]+2.5,
               linetype="solid", color="gray25") +
  geom_segment(y=-0.02, x=Breaks[3]+2.5, xend=Breaks[4],
               linetype="solid", color="gray25") +
  geom_segment(y=-0.0105, x=Breaks[4], xend=Breaks[5],
               linetype="solid", color="gray25") +
  geom_segment(y=-0.01428571, x=Breaks[5], xend=Breaks[6],
               linetype="solid", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
  theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  scale_y_continuous(limits=c(-0.04,0.02),
                     breaks=seq(-0.04,0.02,by=0.01),
                     expand=c(0,0)) +
  scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
  geom_vline(xintercept=Breaks)

Est9 <-
  ggplot(res_df %>% filter(Scenario==9,
                           Estimand %in% c("Overall","Exp. Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_segment(y=-0.01, x=Breaks[3], xend=Breaks[3]+2.5,
               linetype="solid", color="gray25") +
  geom_segment(y=-0.02, x=Breaks[3]+2.5, xend=Breaks[4],
               linetype="solid", color="gray25") +
  geom_segment(y=-0.045, x=Breaks[4], xend=Breaks[5],
               linetype="solid", color="gray25") +
  geom_segment(y=-0.036666667, x=Breaks[5], xend=Breaks[6],
               linetype="solid", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
  theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  scale_y_continuous(limits=c(-0.10,0.0),
                     breaks=seq(-0.10,0.00,by=0.01),
                     expand=c(0,0)) +
  scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
  geom_vline(xintercept=Breaks)

ggsave(filename=paste0(outdir,paste0("Sim_Ests_",outname,".eps")),
       plot=Est1 + guides(color="none", shape="none") + labs(x=NULL, title="A) Scenario 1") +
         Est2 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="B) Scenario 2") +
         Est3 + guides(shape="none") + labs(x=NULL, y=NULL, title="C) Scenario 3") +
         Est4 + guides(color="none", shape="none") + labs(x=NULL, title="D) Scenario 4") +
         Est5 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="E) Scenario 5") +
         Est6 + guides(color="none") + labs(x=NULL, y=NULL, title="F) Scenario 6") +
         Est7 + guides(color="none", shape="none") + labs(title="G) Scenario 7") +
         Est8 + guides(color="none", shape="none") + labs(y=NULL, title="H) Scenario 8") +
         Est9 + guides(color="none", shape="none") + labs(y=NULL, title="I) Scenario 9") +
         plot_layout(nrow=3, ncol=3, byrow=TRUE, guides="collect") +
         theme(legend.position="bottom"),
       width=13, height=12, units="in")
}

CI_Plots <- function(res_df, outname,
                        Labels, Breaks,
                     EstsR1, EstsR2, EstsR3) {

CI1 <-
  ggplot(res_df  %>% filter(Scenario==1,
                            Estimand %in% EstsR1),
         mapping=aes(x=`Estimator Number`, y=`Mean CI Width`,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Mean 95% CI Width") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,0.1),
                     breaks=seq(0,0.1,by=0.02),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

CI2 <-
  ggplot(res_df  %>% filter(Scenario==2,
                            Estimand %in% EstsR1),
         mapping=aes(x=`Estimator Number`, y=`Mean CI Width`,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Mean 95% CI Width") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,0.1),
                     breaks=seq(0,0.1,by=0.02),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

CI3 <-
  ggplot(res_df  %>% filter(Scenario==3,
                            Estimand %in% EstsR1),
         mapping=aes(x=`Estimator Number`, y=`Mean CI Width`,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Mean 95% CI Width") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,0.1),
                     breaks=seq(0,0.1,by=0.02),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

CI4 <-
  ggplot(res_df  %>% filter(Scenario==4,
                            Estimand %in% EstsR2),
         mapping=aes(x=`Estimator Number`, y=`Mean CI Width`,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Mean 95% CI Width") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,0.1),
                     breaks=seq(0,0.1,by=0.02),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

CI5 <-
  ggplot(res_df  %>% filter(Scenario==5,
                            Estimand %in% EstsR2),
         mapping=aes(x=`Estimator Number`, y=`Mean CI Width`,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Mean 95% CI Width") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,0.1),
                     breaks=seq(0,0.1,by=0.02),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

CI6 <-
  ggplot(res_df  %>% filter(Scenario==6,
                            Estimand %in% EstsR2),
         mapping=aes(x=`Estimator Number`, y=`Mean CI Width`,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Mean 95% CI Width") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,0.1),
                     breaks=seq(0,0.1,by=0.02),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

CI7 <-
  ggplot(res_df  %>% filter(Scenario==7,
                            Estimand %in% EstsR3),
         mapping=aes(x=`Estimator Number`, y=`Mean CI Width`,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Mean 95% CI Width") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,0.1),
                     breaks=seq(0,0.1,by=0.02),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

CI8 <-
  ggplot(res_df  %>% filter(Scenario==8,
                            Estimand %in% EstsR3),
         mapping=aes(x=`Estimator Number`, y=`Mean CI Width`,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Mean 95% CI Width") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,0.1),
                     breaks=seq(0,0.1,by=0.02),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

CI9 <-
  ggplot(res_df  %>% filter(Scenario==9,
                            Estimand %in% EstsR3),
         mapping=aes(x=`Estimator Number`, y=`Mean CI Width`,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="Mean 95% CI Width") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(0,0.1),
                     breaks=seq(0,0.1,by=0.02),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

ggsave(filename=paste0(outdir,paste0("Sim_CIs_",outname,".eps")),
       plot=CI1 + guides(color="none", shape="none") + labs(x=NULL, title="A) Scenario 1") +
         CI2 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="B) Scenario 2") +
         CI3 + guides(shape="none") + labs(x=NULL, y=NULL, title="C) Scenario 3") +
         CI4 + guides(color="none", shape="none") + labs(x=NULL, title="D) Scenario 4") +
         CI5 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="E) Scenario 5") +
         CI6 + guides(color="none") + labs(x=NULL, y=NULL, title="F) Scenario 6") +
         CI7 + guides(color="none", shape="none") + labs(title="G) Scenario 7") +
         CI8 + guides(color="none", shape="none") + labs(y=NULL, title="H) Scenario 8") +
         CI9 + guides(color="none", shape="none") + labs(y=NULL, title="I) Scenario 9") +
         plot_layout(nrow=3, ncol=3, byrow=TRUE, guides="collect") +
         theme(legend.position="bottom"),
       width=13, height=12, units="in")
}

Coverage_Plots <- function(res_df, outname,
                           Labels, Breaks,
                           EstsR1, EstsR2, EstsR3) {
Cov1 <-
  ggplot(res_df  %>% filter(Scenario==1,
                            Estimand %in% EstsR1),
         mapping=aes(x=`Estimator Number`, y=100*`CI Coverage`,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=c(95-3*100*sqrt(.95*.05/1000),
                          95+3*100*sqrt(.95*.05/1000)),
             linetype="dotted", color="gray25") +
  geom_hline(yintercept=95, linetype="dashed", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="95% CI Coverage (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(90,100),
                     breaks=seq(90,100,by=2.5),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Cov2 <-
  ggplot(res_df  %>% filter(Scenario==2,
                            Estimand %in% EstsR1),
         mapping=aes(x=`Estimator Number`, y=100*`CI Coverage`,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=c(95-3*100*sqrt(.95*.05/1000),
                          95+3*100*sqrt(.95*.05/1000)),
             linetype="dotted", color="gray25") +
  geom_hline(yintercept=95, linetype="dashed", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="95% CI Coverage (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(90,100),
                     breaks=seq(90,100,by=2.5),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Cov3 <-
  ggplot(res_df  %>% filter(Scenario==3,
                            Estimand %in% EstsR1),
         mapping=aes(x=`Estimator Number`, y=100*`CI Coverage`,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=c(95-3*100*sqrt(.95*.05/1000),
                          95+3*100*sqrt(.95*.05/1000)),
             linetype="dotted", color="gray25") +
  geom_hline(yintercept=95, linetype="dashed", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="95% CI Coverage (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(90,100),
                     breaks=seq(90,100,by=2.5),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Cov4 <-
  ggplot(res_df  %>% filter(Scenario==4,
                            Estimand %in% EstsR2),
         mapping=aes(x=`Estimator Number`, y=100*`CI Coverage`,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=c(95-3*100*sqrt(.95*.05/1000),
                          95+3*100*sqrt(.95*.05/1000)),
             linetype="dotted", color="gray25") +
  geom_hline(yintercept=95, linetype="dashed", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="95% CI Coverage (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(90,100),
                     breaks=seq(90,100,by=2.5),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Cov5 <-
  ggplot(res_df  %>% filter(Scenario==5,
                            Estimand %in% EstsR2),
         mapping=aes(x=`Estimator Number`, y=100*`CI Coverage`,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=c(95-3*100*sqrt(.95*.05/1000),
                          95+3*100*sqrt(.95*.05/1000)),
             linetype="dotted", color="gray25") +
  geom_hline(yintercept=95, linetype="dashed", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="95% CI Coverage (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(90,100),
                     breaks=seq(90,100,by=2.5),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Cov6 <-
  ggplot(res_df  %>% filter(Scenario==6,
                            Estimand %in% EstsR2),
         mapping=aes(x=`Estimator Number`, y=100*`CI Coverage`,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=c(95-3*100*sqrt(.95*.05/1000),
                          95+3*100*sqrt(.95*.05/1000)),
             linetype="dotted", color="gray25") +
  geom_hline(yintercept=95, linetype="dashed", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="95% CI Coverage (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(90,100),
                     breaks=seq(90,100,by=2.5),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Cov7 <-
  ggplot(res_df  %>% filter(Scenario==7,
                            Estimand %in% EstsR3),
         mapping=aes(x=`Estimator Number`, y=100*`CI Coverage`,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=c(95-3*100*sqrt(.95*.05/1000),
                          95+3*100*sqrt(.95*.05/1000)),
             linetype="dotted", color="gray25") +
  geom_hline(yintercept=95, linetype="dashed", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="95% CI Coverage (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(90,100),
                     breaks=seq(90,100,by=2.5),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Cov8 <-
  ggplot(res_df  %>% filter(Scenario==8,
                            Estimand %in% EstsR3),
         mapping=aes(x=`Estimator Number`, y=100*`CI Coverage`,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=c(95-3*100*sqrt(.95*.05/1000),
                          95+3*100*sqrt(.95*.05/1000)),
             linetype="dotted", color="gray25") +
  geom_hline(yintercept=95, linetype="dashed", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="95% CI Coverage (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(90,100),
                     breaks=seq(90,100,by=2.5),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Cov9 <-
  ggplot(res_df  %>% filter(Scenario==9,
                            Estimand %in% EstsR3),
         mapping=aes(x=`Estimator Number`, y=100*`CI Coverage`,
                     color=Type, shape=Assumption)) +
  geom_hline(yintercept=c(95-3*100*sqrt(.95*.05/1000),
                          95+3*100*sqrt(.95*.05/1000)),
             linetype="dotted", color="gray25") +
  geom_hline(yintercept=95, linetype="dashed", color="gray25") +
  geom_point(show.legend=TRUE, size=3) + theme_bw() +
  labs(x="Estimator",y="95% CI Coverage (%)") +
  scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                     expand=expansion(0,0),
                     breaks=Labels$Xval,
                     labels=Labels$Name,
                     minor_breaks=Breaks) +
  geom_vline(xintercept=Breaks) +
  scale_y_continuous(limits=c(90,100),
                     breaks=seq(90,100,by=2.5),
                     expand=c(0,0)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

### Export Plots
ggsave(filename=paste0(outdir,paste0("Sim_Cov_",outname,".eps")),
       plot=Cov1 + guides(color="none", shape="none") + labs(x=NULL, title="A) Scenario 1") +
         Cov2 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="B) Scenario 2") +
         Cov3 + guides(shape="none") + labs(x=NULL, y=NULL, title="C) Scenario 3") +
         Cov4 + guides(color="none", shape="none") + labs(x=NULL, title="D) Scenario 4") +
         Cov5 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="E) Scenario 5") +
         Cov6 + guides(color="none") + labs(x=NULL, y=NULL, title="F) Scenario 6") +
         Cov7 + guides(color="none", shape="none") + labs(title="G) Scenario 7") +
         Cov8 + guides(color="none", shape="none") + labs(y=NULL, title="H) Scenario 8") +
         Cov9 + guides(color="none", shape="none") + labs(y=NULL, title="I) Scenario 9") +
         plot_layout(nrow=3, ncol=3, byrow=TRUE, guides="collect") +
         theme(legend.position="bottom"),
       width=13, height=12, units="in")
}

Target_Est_Plots <- function(res_df, outname,
                              Labels, Breaks) {
  ### Plots of Estimates and SDs:
  Est1 <-
    ggplot(res_df %>% filter(Scenario==1),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_hline(yintercept=0, linetype="solid", color="gray25") +
    geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
    theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                       expand=expansion(0,0),
                       breaks=Labels$Xval,
                       labels=Labels$Name,
                       minor_breaks=Breaks) +
    scale_y_continuous(limits=c(-0.02,0.02),
                       breaks=seq(-0.02,0.02,by=0.005),
                       expand=c(0,0)) +
    scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
    geom_vline(xintercept=Breaks)

  Est2 <-
    ggplot(res_df %>% filter(Scenario==2),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_hline(yintercept=-0.02, linetype="solid", color="gray25") +
    geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
    theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                       expand=expansion(0,0),
                       breaks=Labels$Xval,
                       labels=Labels$Name,
                       minor_breaks=Breaks) +
    scale_y_continuous(limits=c(-0.04,0.00),
                       breaks=seq(-0.04,0.00,by=0.005),
                       expand=c(0,0)) +
    scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
    geom_vline(xintercept=Breaks)

  Est3 <-
    ggplot(res_df %>% filter(Scenario==3),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_hline(yintercept=-0.04, linetype="solid", color="gray25") +
    geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
    theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                       expand=expansion(0,0),
                       breaks=Labels$Xval,
                       labels=Labels$Name,
                       minor_breaks=Breaks) +
    scale_y_continuous(limits=c(-0.06, -0.02),
                       breaks=seq(-0.06, -0.02,by=0.005),
                       expand=c(0,0)) +
    scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
    geom_vline(xintercept=Breaks)

  Est4 <-
    ggplot(res_df %>% filter(Scenario==4,
                             Estimand=="Time 3"),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_segment(y=-0.05, x=Breaks[1], xend=Breaks[2],
                 linetype="solid", color="gray25") +
    geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
    theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                       expand=expansion(0,0),
                       breaks=Labels$Xval,
                       labels=Labels$Name,
                       minor_breaks=Breaks) +
    scale_y_continuous(limits=c(-0.08,0),
                       breaks=seq(-0.08,0,by=0.01),
                       expand=c(0,0)) +
    scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
    geom_vline(xintercept=Breaks)

  Est5 <-
    ggplot(res_df %>% filter(Scenario==5,
                             Estimand=="Time 3"),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_segment(y=-0.06, x=Breaks[1], xend=Breaks[2],
                 linetype="solid", color="gray25") +
    geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
    theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                       expand=expansion(0,0),
                       breaks=Labels$Xval,
                       labels=Labels$Name,
                       minor_breaks=Breaks) +
    scale_y_continuous(limits=c(-0.08,0),
                       breaks=seq(-0.08,0,by=0.01),
                       expand=c(0,0)) +
    scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
    geom_vline(xintercept=Breaks)

  Est6 <-
    ggplot(res_df %>% filter(Scenario==6,
                             Estimand=="Time 3"),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_segment(y=-0.03, x=Breaks[1], xend=Breaks[2],
                 linetype="solid", color="gray25") +
    geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
    theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                       expand=expansion(0,0),
                       breaks=Labels$Xval,
                       labels=Labels$Name,
                       minor_breaks=Breaks) +
    scale_y_continuous(limits=c(-0.08,0),
                       breaks=seq(-0.08,0,by=0.01),
                       expand=c(0,0)) +
    scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
    geom_vline(xintercept=Breaks)

  Est7 <-
    ggplot(res_df %>% filter(Scenario==7,
                             Estimand %in% c("Exp. Pd. 1","Exp. Pd. 2")),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_segment(y=-0.015, x=Breaks[2], xend=Breaks[3],
                 linetype="solid", color="gray25") +
    geom_segment(y=-0.010, x=Breaks[3], xend=Breaks[4],
                 linetype="solid", color="gray25") +
    geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
    theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                       expand=expansion(0,0),
                       breaks=Labels$Xval,
                       labels=Labels$Name,
                       minor_breaks=Breaks) +
    scale_y_continuous(limits=c(-0.04,0.02),
                       breaks=seq(-0.04,0.02,by=0.01),
                       expand=c(0,0)) +
    scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
    geom_vline(xintercept=Breaks)

  Est8 <-
    ggplot(res_df %>% filter(Scenario==8,
                             Estimand %in% c("Exp. Pd. 1","Exp. Pd. 2")),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_segment(y=0, x=Breaks[2], xend=Breaks[3],
                 linetype="solid", color="gray25") +
    geom_segment(y=0, x=Breaks[3], xend=Breaks[4],
                 linetype="solid", color="gray25") +
    geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
    theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                       expand=expansion(0,0),
                       breaks=Labels$Xval,
                       labels=Labels$Name,
                       minor_breaks=Breaks) +
    scale_y_continuous(limits=c(-0.04,0.02),
                       breaks=seq(-0.04,0.02,by=0.01),
                       expand=c(0,0)) +
    scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
    geom_vline(xintercept=Breaks)

  Est9 <-
    ggplot(res_df %>% filter(Scenario==9,
                             Estimand %in% c("Exp. Pd. 1","Exp. Pd. 2")),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_segment(y=-0.05, x=Breaks[2], xend=Breaks[3],
                 linetype="solid", color="gray25") +
    geom_segment(y=-0.07, x=Breaks[3], xend=Breaks[4],
                 linetype="solid", color="gray25") +
    geom_point(show.legend=TRUE, size=3) + geom_errorbar(linewidth=1.2) +
    theme_bw() + labs(x="Estimand", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(Breaks[1], Breaks[length(Breaks)]),
                       expand=expansion(0,0),
                       breaks=Labels$Xval,
                       labels=Labels$Name,
                       minor_breaks=Breaks) +
    scale_y_continuous(limits=c(-0.10,-0.02),
                       breaks=seq(-0.10,-0.02,by=0.01),
                       expand=c(0,0)) +
    scale_color_manual(name="Estimator Type", drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(name="Assumption", drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2")) +
    geom_vline(xintercept=Breaks)

  ggsave(filename=paste0(outdir,paste0("Sim_Ests_",outname,".eps")),
         plot=Est1 + guides(color="none", shape="none") + labs(x=NULL, title="A) Scenario 1") +
           Est2 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="B) Scenario 2") +
           Est3 + guides(shape="none") + labs(x=NULL, y=NULL, title="C) Scenario 3") +
           Est4 + guides(color="none", shape="none") + labs(x=NULL, title="D) Scenario 4") +
           Est5 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="E) Scenario 5") +
           Est6 + guides(color="none") + labs(x=NULL, y=NULL, title="F) Scenario 6") +
           Est7 + guides(color="none", shape="none") + labs(title="G) Scenario 7") +
           Est8 + guides(color="none", shape="none") + labs(y=NULL, title="H) Scenario 8") +
           Est9 + guides(color="none", shape="none") + labs(y=NULL, title="I) Scenario 9") +
           plot_layout(nrow=3, ncol=3, byrow=TRUE, guides="collect") +
           theme(legend.position="bottom"),
         width=13, height=12, units="in")
}

## Get manuscript plots for power, CI width, and CI coverage:

Power_Plots(Overall, outname="Overall_Ind",
            Labels=Overall_Labels, Breaks=Overall_Breaks,
            EstsR1=c("Overall", "Time Avg.", "Exp. Avg.","Group Avg.","ATT"),
            EstsR2=c("Time Avg.","Group Avg.","ATT"),
            EstsR3=c("Exp. Avg.","Group Avg.","ATT"))

CI_Plots(Overall, outname="Overall_Ind",
            Labels=Overall_Labels, Breaks=Overall_Breaks,
         EstsR1=c("Overall", "Time Avg.", "Exp. Avg.","Group Avg.","ATT"),
         EstsR2=c("Overall", "Time Avg.","Group Avg.","ATT"),
         EstsR3=c("Overall", "Exp. Avg.","Group Avg.","ATT"))

Coverage_Plots(Overall, outname="Overall_Ind",
            Labels=Overall_Labels, Breaks=Overall_Breaks,
            EstsR1=c("Overall", "Time Avg.", "Exp. Avg.","Group Avg.","ATT"),
            EstsR2=c("Time Avg.","Group Avg.","ATT"),
            EstsR3=c("Exp. Avg.","Group Avg.","ATT"))

Power_Plots(Target, outname="Target_Ind",
            Labels=Target_Labels, Breaks=Target_Breaks,
            EstsR1=c("Time 3","Exp. Pd. 2","Exp. Pd. 1"),
            EstsR2=c("Time 3"),
            EstsR3=c("Exp. Pd. 2","Exp. Pd. 1"))

CI_Plots(Target, outname="Target_Ind",
         Labels=Target_Labels, Breaks=Target_Breaks,
         EstsR1=c("Time 3","Exp. Pd. 2","Exp. Pd. 1"),
         EstsR2=c("Time 3"),
         EstsR3=c("Exp. Pd. 2","Exp. Pd. 1"))

Coverage_Plots(Target, outname="Target_Ind",
               Labels=Target_Labels, Breaks=Target_Breaks,
               EstsR1=c("Time 3","Exp. Pd. 2","Exp. Pd. 1"),
               EstsR2=c("Time 3"),
               EstsR3=c("Exp. Pd. 2","Exp. Pd. 1"))

Power_Plots(Overall_CS, outname="Overall_CS",
            Labels=Overall_Labels, Breaks=Overall_Breaks,
            EstsR1=c("Overall", "Time Avg.", "Exp. Avg.","Group Avg.","ATT"),
            EstsR2=c("Time Avg.","Group Avg.","ATT"),
            EstsR3=c("Exp. Avg.","Group Avg.","ATT"))

CI_Plots(Overall_CS, outname="Overall_CS",
         Labels=Overall_Labels, Breaks=Overall_Breaks,
         EstsR1=c("Overall", "Time Avg.", "Exp. Avg.","Group Avg.","ATT"),
         EstsR2=c("Overall", "Time Avg.","Group Avg.","ATT"),
         EstsR3=c("Overall", "Exp. Avg.","Group Avg.","ATT"))

Coverage_Plots(Overall_CS, outname="Overall_CS",
               Labels=Overall_Labels, Breaks=Overall_Breaks,
               EstsR1=c("Overall", "Time Avg.", "Exp. Avg.","Group Avg.","ATT"),
               EstsR2=c("Time Avg.","Group Avg.","ATT"),
               EstsR3=c("Exp. Avg.","Group Avg.","ATT"))

Power_Plots(Target_CS, outname="Target_CS",
            Labels=Target_Labels, Breaks=Target_Breaks,
            EstsR1=c("Time 3","Exp. Pd. 2","Exp. Pd. 1"),
            EstsR2=c("Time 3"),
            EstsR3=c("Exp. Pd. 2","Exp. Pd. 1"))

CI_Plots(Target_CS, outname="Target_CS",
         Labels=Target_Labels, Breaks=Target_Breaks,
         EstsR1=c("Time 3","Exp. Pd. 2","Exp. Pd. 1"),
         EstsR2=c("Time 3"),
         EstsR3=c("Exp. Pd. 2","Exp. Pd. 1"))

Coverage_Plots(Target_CS, outname="Target_CS",
               Labels=Target_Labels, Breaks=Target_Breaks,
               EstsR1=c("Time 3","Exp. Pd. 2","Exp. Pd. 1"),
               EstsR2=c("Time 3"),
               EstsR3=c("Exp. Pd. 2","Exp. Pd. 1"))

## Get manuscript plots for estimates for Overall set:

Overall_Est_Plots(Overall, outname="Overall_Ind",
                  Labels=Overall_Labels, Breaks=Overall_Breaks)

Overall_Est_Plots(Overall_CS, outname="Overall_CS",
                  Labels=Overall_Labels, Breaks=Overall_Breaks)

## Get manuscript plots for estimates for Target set:

Target_Est_Plots(Target, outname="Target_Ind",
                  Labels=Target_Labels, Breaks=Target_Breaks)

Target_Est_Plots(Target_CS, outname="Target_CS",
                  Labels=Target_Labels, Breaks=Target_Breaks)
