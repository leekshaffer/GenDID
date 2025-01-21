#######################################
###### File: Sim_Plots.R ##############
###### Lee Kennedy-Shaffer ############
###### Created 2024/10/03 #############
#######################################

require(tidyverse)
require(patchwork)

Colors <- c("#1b9e77","#d95f02","#7570b3","#e7298a")
Shapes <- c(15:18,8)

outdir <- "figs/"

## Load simulation results:
load(file="int/Full_Sim_Res.Rda")

## Select Desired Estimators:
OverallSet <- tibble(Estimator=c("A5_Ind_","Comp_W_TW","Comp_CPI","Comp_W_CO.W_CO3",
                                 "A4_Ind_AvgEx8","CPI.T_AvgExLast","Comp_W_CS.W_calendar","A2_Ind_T.Avg","CPI.DT_TAvg",
                                 "A3_Ind_Avg","CPI.D_Avg","A3_Ind_AvgEx7","Comp_W_CS.W_dynamic","A2_Ind_D.Avg","CPI.DT_DAvg",
                                 "A2_Ind_Group","Comp_W_CS.W_group","Comp_CLWP","Comp_CLWPA",
                                 "A2_Ind_AvgExT8","Comp_W_CS.W_simple","Comp_W_SA.W_ATT","CPI.DT_AvgEx8")) %>%
  mutate(`Estimator Number`=row_number(),
         Type=c("GD","SA","ME","SA",
                "GD","ME","SA","GD","ME",
                "GD","ME","GD","SA","GD","ME",
                "GD","SA","CL","CL",
                "GD","SA","SA","ME"),
         Assumption=c("S5","S5","S5","S5",
                      "S4","S4","S4","S2","S2","S4","S4",
                      "S3","S3","S3","S3","S2","S2",
                      "S2","S2",
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
OverallSet_333 <- tibble(Estimator=c("A5_333_","Comp_W_TW","Comp_CPI","Comp_W_CO.W_CO3",
                                 "A4_333_AvgEx8","CPI.T_AvgExLast","Comp_W_CS.W_calendar","A2_333_T.Avg","CPI.DT_TAvg",
                                 "A3_333_Avg","CPI.D_Avg","A3_333_AvgEx7","Comp_W_CS.W_dynamic","A2_333_D.Avg","CPI.DT_DAvg",
                                 "A2_333_Group","Comp_W_CS.W_group","Comp_CLWP","Comp_CLWPA",
                                 "A2_333_AvgExT8","Comp_W_CS.W_simple","Comp_W_SA.W_ATT","CPI.DT_AvgEx8")) %>%
  mutate(`Estimator Number`=row_number(),
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
TargetsSet <- tibble(Estimator=c("A4_Ind_T.3","CPI.T_3","A2_Ind_T.3","CPI.DT_T3",
                                 "A3_Ind_D.2","CPI.D_2","A2_Ind_D.2","CPI.DT_D2",
                                 "A3_Ind_D.1","CPI.D_1","Comp_W_CO.W_CO2",
                                 "A2_Ind_D.1","CPI.DT_D1","Comp_W_CH.W_M","Comp_W_CO.W_CO1"
                                 )) %>%
  mutate(`Estimator Number`=row_number(),
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
TargetsSet_333 <- tibble(Estimator=c("A4_333_T.3","CPI.T_3","A2_333_T.3","CPI.DT_T3",
                                 "A3_333_D.2","CPI.D_2","A2_333_D.2","CPI.DT_D2",
                                 "A3_333_D.1","CPI.D_1","Comp_W_CO.W_CO2",
                                 "A2_333_D.1","CPI.DT_D1","Comp_W_CH.W_M","Comp_W_CO.W_CO1"
)) %>%
  mutate(`Estimator Number`=row_number(),
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
Overall <- Full_Sim_Res %>% dplyr::select(all_of(c("SimNo","Result",OverallSet$Estimator))) %>%
  pivot_longer(cols=-c("SimNo","Result"),names_to="Estimator", values_to="Value") %>%
  pivot_wider(id_cols=c("SimNo","Estimator"), names_from="Result", values_from="Value") %>%
  mutate(Lower=`Mean Estimate`-`SD Estimate`, Upper=`Mean Estimate`+`SD Estimate`) %>%
  left_join(OverallSet, by="Estimator")
Overall_333 <- Full_Sim_Res %>% dplyr::select(all_of(c("SimNo","Result",OverallSet_333$Estimator))) %>%
  pivot_longer(cols=-c("SimNo","Result"),names_to="Estimator", values_to="Value") %>%
  pivot_wider(id_cols=c("SimNo","Estimator"), names_from="Result", values_from="Value") %>%
  mutate(Lower=`Mean Estimate`-`SD Estimate`, Upper=`Mean Estimate`+`SD Estimate`) %>%
  left_join(OverallSet_333, by="Estimator")
Targets <- Full_Sim_Res %>% dplyr::select(all_of(c("SimNo","Result",TargetsSet$Estimator))) %>%
  pivot_longer(cols=-c("SimNo","Result"),names_to="Estimator", values_to="Value") %>%
  pivot_wider(id_cols=c("SimNo","Estimator"), names_from="Result", values_from="Value") %>%
  mutate(Lower=`Mean Estimate`-`SD Estimate`, Upper=`Mean Estimate`+`SD Estimate`) %>%
  left_join(TargetsSet, by="Estimator")
Targets_333 <- Full_Sim_Res %>% dplyr::select(all_of(c("SimNo","Result",TargetsSet_333$Estimator))) %>%
  pivot_longer(cols=-c("SimNo","Result"),names_to="Estimator", values_to="Value") %>%
  pivot_wider(id_cols=c("SimNo","Estimator"), names_from="Result", values_from="Value") %>%
  mutate(Lower=`Mean Estimate`-`SD Estimate`, Upper=`Mean Estimate`+`SD Estimate`) %>%
  left_join(TargetsSet_333, by="Estimator")

### Plots for full set:
Overall_Plots <- function(res_df, outname,
                          BreakVec, MinorVec) {
Power1 <-
  ggplot(res_df  %>% filter(SimNo==1),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=2) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0),
                     breaks=BreakVec, minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  geom_hline(yintercept=5, linetype="dotted", color="gray25") +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power2 <-
  ggplot(res_df  %>% filter(SimNo==2),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=2) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0),
                     breaks=BreakVec, minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power3 <-
  ggplot(res_df  %>% filter(SimNo==3),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=2) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0),
                     breaks=BreakVec, minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power4 <-
  ggplot(res_df  %>% filter(SimNo==4, #Assumption %in% c("S4","S5")),
                             Estimand %in% c("Time Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=2) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0),
                     breaks=BreakVec, minor_breaks=MinorVec) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power5 <-
  ggplot(res_df  %>% filter(SimNo==5, #Assumption %in% c("S4","S5")),
                             Estimand %in% c("Time Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=2) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0),
                     breaks=BreakVec, minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power6 <-
  ggplot(res_df  %>% filter(SimNo==6, #Assumption %in% c("S4","S5")),
                             Estimand %in% c("Time Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=2) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power7 <-
  ggplot(res_df  %>% filter(SimNo==7, #Assumption %in% c("S3")),
                             Estimand %in% c("Exp. Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=2) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power8 <-
  ggplot(res_df  %>% filter(SimNo==8, #Assumption %in% c("S3")),
                             Estimand %in% c("Exp. Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=2) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Power9 <-
  ggplot(res_df  %>% filter(SimNo==9, #Assumption %in% c("S3")),
                             Estimand %in% c("Exp. Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=Power*100,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE, size=2) + theme_bw() +
  labs(x="Estimator",y="Empirical Power (%)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

### Plots of Estimates and SDs:
Est1 <-
  ggplot(res_df %>% filter(SimNo==1),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(-0.10,0.015),
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_hline(yintercept=0, linetype="dotted", color="gray25") +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Est2 <-
  ggplot(res_df %>% filter(SimNo==2),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(-0.10,0.015),
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_hline(yintercept=-0.02, linetype="dotted", color="gray25") +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Est3 <-
  ggplot(res_df %>% filter(SimNo==3),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(-0.10,0.015),
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_hline(yintercept=-0.04, linetype="dotted", color="gray25") +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Est4 <-
  ggplot(res_df %>% filter(SimNo==4,
                            # Assumption %in% c("S2","S4","S5")),
                            Estimand %in% c("Overall","Time Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(-0.10,0.015),
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.02, x=0.5, xend=9.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=0.005, x=15.5, xend=17.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.02, x=17.5, xend=19.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.0033333, x=19.5, xend=23.5,
               linetype="dotted", color="gray25") +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Est5 <-
  ggplot(res_df %>% filter(SimNo==5,
                            # Assumption %in% c("S2","S4","S5")),
                            Estimand %in% c("Overall","Time Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(-0.10,0.015),
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.02, x=0.5, xend=9.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=0.005694444, x=15.5, xend=17.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.02, x=17.5, xend=19.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.001904762, x=19.5, xend=23.5,
               linetype="dotted", color="gray25") +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Est6 <-
  ggplot(res_df %>% filter(SimNo==6,
                            # Assumption %in% c("S2","S4","S5")),
                            Estimand %in% c("Overall","Time Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(-0.10,0.015),
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.02, x=0.5, xend=9.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.0105, x=15.5, xend=17.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.02, x=17.5, xend=19.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.01428571, x=19.5, xend=23.5,
               linetype="dotted", color="gray25") +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Est7 <-
  ggplot(res_df %>% filter(SimNo==7,
                            # Assumption %in% c("S2","S3","S5")),
                            Estimand %in% c("Overall","Exp. Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(-0.10,0.015),
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.025, x=9.5, xend=11.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.0225, x=11.5, xend=15.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-.01625, x=15.5, xend=19.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.01833333, x=19.5, xend=23.5,
               linetype="dotted", color="gray25") +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Est8 <-
  ggplot(res_df %>% filter(SimNo==8,
                            # Assumption %in% c("S2","S3","S5")),
                            Estimand %in% c("Overall","Exp. Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(-0.10,0.015),
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.0214, x=9.5, xend=11.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.02, x=11.5, xend=15.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.0105, x=15.5, xend=19.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.01428571, x=19.5, xend=23.5,
               linetype="dotted", color="gray25") +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

Est9 <-
  ggplot(res_df %>% filter(SimNo==9,
                            # Assumption %in% c("S2","S3","S5")),
                            Estimand %in% c("Overall","Exp. Avg.","Group Avg.","ATT")),
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) +
  geom_point(show.legend=TRUE) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_continuous(limits=c(0,24), expand=expansion(0,0), breaks=BreakVec,
                     minor_breaks=MinorVec) +
  scale_y_continuous(limits=c(-0.10,0.015),
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.01, x=9.5, xend=11.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.02, x=11.5, xend=15.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.045, x=15.5, xend=19.5,
               linetype="dotted", color="gray25") +
  geom_segment(y=-0.036666667, x=19.5, xend=23.5,
               linetype="dotted", color="gray25") +
  scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
  scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

### Export Plots
ggsave(filename=paste0(outdir,paste0("Sim_Power_",outname,".eps")),
       plot=Power1 + guides(color="none", shape="none") + labs(x=NULL, title="A) Scenario 1") +
         Power2 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="B) Scenario 2") +
         Power3 + guides(shape="none") +labs(x=NULL, y=NULL, title="C) Scenario 3") +
         Power4 + guides(color="none", shape="none") + labs(x=NULL, title="D) Scenario 4") +
         Power5 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="E) Scenario 5") +
         Power6 + guides(color="none") + labs(x=NULL, y=NULL, title="F) Scenario 6") +
         Power7 + guides(color="none", shape="none") + labs(title="G) Scenario 7") +
         Power8 + guides(color="none", shape="none") + labs(y=NULL, title="H) Scenario 8") +
         Power9 + guides(color="none", shape="none") + labs(y=NULL, title="I) Scenario 9") +
         plot_layout(nrow=3, ncol=3, byrow=TRUE),
       width=8, height=7, units="in")
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
         plot_layout(nrow=3, ncol=3, byrow=TRUE),
       width=8, height=7, units="in")
}

Overall_Plots(Overall, outname="Overall_Ind",
              BreakVec=seq(3,21,by=3), MinorVec=NULL)
Overall_Plots(Overall_333, outname="Overall_333",
              BreakVec=seq(3,21,by=3), MinorVec=NULL)

### Plots for Targeted period-specific effects:
Target_Plots <- function(res_df, outname,
                         BreakVec, MinorVec) {
  Power1 <-
    ggplot(res_df  %>% filter(SimNo==1),
           mapping=aes(x=`Estimator Number`, y=Power*100,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE, size=2) + theme_bw() +
    labs(x="Estimator",y="Empirical Power (%)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(0,100),
                       breaks=seq(0,100,by=20),
                       expand=c(0,2)) +
    geom_hline(yintercept=5, linetype="dotted", color="gray25") +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Power2 <-
    ggplot(res_df  %>% filter(SimNo==2),
           mapping=aes(x=`Estimator Number`, y=Power*100,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE, size=2) + theme_bw() +
    labs(x="Estimator",y="Empirical Power (%)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(0,100),
                       breaks=seq(0,100,by=20),
                       expand=c(0,2)) +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Power3 <-
    ggplot(res_df  %>% filter(SimNo==3),
           mapping=aes(x=`Estimator Number`, y=Power*100,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE, size=2) + theme_bw() +
    labs(x="Estimator",y="Empirical Power (%)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(0,100),
                       breaks=seq(0,100,by=20),
                       expand=c(0,2)) +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Power4 <-
    ggplot(res_df  %>% filter(SimNo==4,
                              Estimand=="Time 3"),
           mapping=aes(x=`Estimator Number`, y=Power*100,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE, size=2) + theme_bw() +
    labs(x="Estimator",y="Empirical Power (%)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(0,100),
                       breaks=seq(0,100,by=20),
                       expand=c(0,2)) +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Power5 <-
    ggplot(res_df  %>% filter(SimNo==5,
                              Estimand=="Time 3"),
           mapping=aes(x=`Estimator Number`, y=Power*100,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE, size=2) + theme_bw() +
    labs(x="Estimator",y="Empirical Power (%)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(0,100),
                       breaks=seq(0,100,by=20),
                       expand=c(0,2)) +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Power6 <-
    ggplot(res_df  %>% filter(SimNo==6,
                              Estimand=="Time 3"),
           mapping=aes(x=`Estimator Number`, y=Power*100,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE, size=2) + theme_bw() +
    labs(x="Estimator",y="Empirical Power (%)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(0,100),
                       breaks=seq(0,100,by=20),
                       expand=c(0,2)) +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Power7 <- ggplot(res_df  %>% filter(SimNo==7,
                                      Estimand %in% c("Exp. Pd. 1","Exp. Pd. 2")),
                   mapping=aes(x=`Estimator Number`, y=Power*100,
                               color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE, size=2) + theme_bw() +
    labs(x="Estimator",y="Empirical Power (%)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(0,100),
                       breaks=seq(0,100,by=20),
                       expand=c(0,2)) +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Power8 <- ggplot(res_df  %>% filter(SimNo==8,
                                      Estimand %in% c("Exp. Pd. 1","Exp. Pd. 2")),
                   mapping=aes(x=`Estimator Number`, y=Power*100,
                               color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE, size=2) + theme_bw() +
    labs(x="Estimator",y="Empirical Power (%)") +
    geom_hline(yintercept=5, linetype="dotted", color="gray25") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(0,100),
                       breaks=seq(0,100,by=20),
                       expand=c(0,2)) +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Power9 <- ggplot(res_df  %>% filter(SimNo==9,
                                      Estimand %in% c("Exp. Pd. 1","Exp. Pd. 2")),
                   mapping=aes(x=`Estimator Number`, y=Power*100,
                               color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE, size=2) + theme_bw() +
    labs(x="Estimator",y="Empirical Power (%)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(0,100),
                       breaks=seq(0,100,by=20),
                       expand=c(0,2)) +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Est1 <-
    ggplot(res_df %>% filter(SimNo==1),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE) + geom_errorbar() +
    theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(-0.10,0.015),
                       breaks=seq(-0.1,0.01,by=0.02),
                       expand=c(0,0)) +
    geom_hline(yintercept=0, linetype="dotted", color="gray25") +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Est2 <-
    ggplot(res_df %>% filter(SimNo==2),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE) + geom_errorbar() +
    theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(-0.10,0.015),
                       breaks=seq(-0.1,0.01,by=0.02),
                       expand=c(0,0)) +
    geom_hline(yintercept=-0.02, linetype="dotted", color="gray25") +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Est3 <-
    ggplot(res_df %>% filter(SimNo==3),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE) + geom_errorbar() +
    theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(-0.10,0.015),
                       breaks=seq(-0.1,0.01,by=0.02),
                       expand=c(0,0)) +
    geom_hline(yintercept=-0.04, linetype="dotted", color="gray25") +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Est4 <-
    ggplot(res_df %>% filter(SimNo==4,
                             # Assumption %in% c("S2","S4","S5")),
                             Estimand=="Time 3"),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE) + geom_errorbar() +
    theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(-0.10,0.015),
                       breaks=seq(-0.1,0.01,by=0.02),
                       expand=c(0,0)) +
    geom_segment(y=-0.05, x=0.5, xend=4.5,
                 linetype="dotted", color="gray25") +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Est5 <-
    ggplot(res_df %>% filter(SimNo==5,
                             # Assumption %in% c("S2","S4","S5")),
                             Estimand=="Time 3"),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE) + geom_errorbar() +
    theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(-0.10,0.015),
                       breaks=seq(-0.1,0.01,by=0.02),
                       expand=c(0,0)) +
    geom_segment(y=-0.06, x=0.5, xend=4.5,
                 linetype="dotted", color="gray25") +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Est6 <-
    ggplot(res_df %>% filter(SimNo==6,
                             # Assumption %in% c("S2","S4","S5")),
                             Estimand=="Time 3"),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE) + geom_errorbar() +
    theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(-0.10,0.015),
                       breaks=seq(-0.1,0.01,by=0.02),
                       expand=c(0,0)) +
    geom_segment(y=-0.03, x=0.5, xend=4.5,
                 linetype="dotted", color="gray25") +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Est7 <-
    ggplot(res_df %>% filter(SimNo==7,
                             # Assumption %in% c("S2","S4","S5")),
                             Estimand %in% c("Exp. Pd. 1","Exp. Pd. 2")),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE) + geom_errorbar() +
    theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(-0.10,0.015),
                       breaks=seq(-0.1,0.01,by=0.02),
                       expand=c(0,0)) +
    geom_segment(y=-0.015, x=4.5, xend=8.5,
                 linetype="dotted", color="gray25") +
    geom_segment(y=-0.01, x=8.5, xend=15.5,
                 linetype="dotted", color="gray25")  +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Est8 <-
    ggplot(res_df %>% filter(SimNo==8,
                             # Assumption %in% c("S2","S4","S5")),
                             Estimand %in% c("Exp. Pd. 1","Exp. Pd. 2")),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE) + geom_errorbar() +
    theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(-0.10,0.015),
                       breaks=seq(-0.1,0.01,by=0.02),
                       expand=c(0,0)) +
    geom_segment(y=0, x=4.5, xend=8.5,
                 linetype="dotted", color="gray25") +
    geom_segment(y=0, x=8.5, xend=15.5,
                 linetype="dotted", color="gray25") +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  Est9 <-
    ggplot(res_df %>% filter(SimNo==9,
                             # Assumption %in% c("S2","S4","S5")),
                             Estimand %in% c("Exp. Pd. 1","Exp. Pd. 2")),
           mapping=aes(x=`Estimator Number`, y=`Mean Estimate`,
                       ymin=Lower, ymax=Upper,
                       color=Type, shape=Assumption)) +
    geom_point(show.legend=TRUE) + geom_errorbar() +
    theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
    scale_x_continuous(limits=c(0,16), expand=expansion(0,0), breaks=BreakVec,
                       minor_breaks=MinorVec) +
    scale_y_continuous(limits=c(-0.10,0.015),
                       breaks=seq(-0.1,0.01,by=0.02),
                       expand=c(0,0)) +
    geom_segment(y=-0.05, x=4.5, xend=8.5,
                 linetype="dotted", color="gray25") +
    geom_segment(y=-0.07, x=8.5, xend=15.5,
                 linetype="dotted", color="gray25")  +
    scale_color_manual(drop=FALSE, limits=c("GD","ME","SA","CL"), values=Colors, breaks=c("GD","ME","SA","CL")) +
    scale_shape_manual(drop=FALSE, limits=c("S5","S4","S3","S2"), values=Shapes, breaks=c("S5","S4","S3","S2"))

  ggsave(filename=paste0(outdir,paste0("Sim_Power_",outname,".eps")),
         plot=Power1 + guides(color="none", shape="none") + labs(x=NULL, title="A) Scenario 1") +
           Power2 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="B) Scenario 2") +
           Power3 + guides(shape="none") +labs(x=NULL, y=NULL, title="C) Scenario 3") +
           Power4 + guides(color="none", shape="none") + labs(x=NULL, title="D) Scenario 4") +
           Power5 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="E) Scenario 5") +
           Power6 + guides(color="none") + labs(x=NULL, y=NULL, title="F) Scenario 6") +
           Power7 + guides(color="none", shape="none") + labs(title="G) Scenario 7") +
           Power8 + guides(color="none", shape="none") + labs(y=NULL, title="H) Scenario 8") +
           Power9 + guides(color="none", shape="none") + labs(y=NULL, title="I) Scenario 9") +
           plot_layout(nrow=3, ncol=3, byrow=TRUE),
         width=8, height=7, units="in")
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
           plot_layout(nrow=3, ncol=3, byrow=TRUE),
         width=8, height=7, units="in")
}

Target_Plots(Targets, outname="Target_Ind",
             BreakVec=seq(3,15,by=3), MinorVec=NULL)
Target_Plots(Targets_333, outname="Target_333",
             BreakVec=seq(3,15,by=3), MinorVec=NULL)
