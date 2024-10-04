#######################################
###### File: Sim_Plots.R ##############
###### Lee Kennedy-Shaffer ############
###### Created 2024/10/03 #############
#######################################

require(tidyverse)
require(patchwork)

outdir <- "figs/"

## Load simulation results from Runs 1--2:
load(file="sim_res/Full_Sim_Res_12.Rda")

## Select Desired Estimators:
OverallSet <- tibble(Estimator=c("A5_Ind_","Comp_W_TW","Comp_CPI","Comp_W_CO.W_CO3",
                                 "A4_Ind_AvgEx8","CPI.T_AvgExLast","Comp_W_CS.W_calendar"
                                 "A3_Ind_Avg","CPI.D_Avg","A3_Ind_AvgEx7","Comp_W_CS.W_dynamic",
                                 "A2_Ind_Group","Comp_W_CS.W_group",
                                 "A2_Ind_AvgEx7","Comp_W_CS.W_simple","Comp_W_SA.W_ATT")) %>%
  ### When update sims, change A2_Ind_AvgEx7 to A2_Ind_AvgExT8 and add A2_Ind_D.Avg in row 3 and A2_Ind_T.Avg in row 2
  mutate(`Estimator Number`=row_number(),
         Type=c("GD","SA","CPI","SA","GD","CPI","SA","GD","CPI","GD","SA",
                "GD","SA","GD","SA","SA"),
         Assumption=c("S5","S5","S5","S5","S4","S4","S4","S3","S3","S3","S3",
                      "S2","S2","S2","S2","S2"),
         Name=c("GD_A5","TWFE","CPI_A5","CO3","GD_A4","CPI_A4","CS_A4",
                "GD_A3","CPI_A3","GD_A3_ExLast","CS_A3",
                "GD_A2_group","CS_group","GD_A2_ATT","CS_ATT","SA_ATT"))
OverallSet_333 <- tibble(Estimator=c("A5_333_","Comp_W_TW","Comp_CPI","Comp_W_CO.W_CO3",
                                     "A4_333_AvgEx8","CPI.T_AvgExLast","Comp_W_CS.W_calendar",
                                     "A3_333_Avg","CPI.D_Avg","A3_333_AvgEx7","Comp_W_CS.W_dynamic",
                                     "A2_333_Group","Comp_W_CS.W_group",
                                     "A2_333_AvgEx7","Comp_W_CS.W_simple","Comp_W_SA.W_ATT")) %>%
  mutate(`Estimator Number`=row_number(),
         Type=c("GD","SA","CPI","SA","GD","CPI","SA","GD","CPI","GD","SA",
                "GD","SA","GD","SA","SA"),
         Assumption=c("S5","S5","S5","S5","S4","S4","S4","S3","S3","S3","S3",
                      "S2","S2","S2","S2","S2"),
         Name=c("GD_A5","TWFE","CPI_A5","CO3","GD_A4","CPI_A4","CS_A4",
                "GD_A3","CPI_A3","GD_A3_ExLast","CS_A3",
                "GD_A2_group","CS_group","GD_A2_ATT","CS_ATT","SA_ATT"))
TargetsSet <- tibble(Estimator=c("A4_Ind_T.3","CPI.T_3","A2_Ind_T.3",
                                 "A3_Ind_D.2","CPI.D_2","A2_Ind_D.2",
                                 "A2_Ind_D.1","Comp_W_CH.W_M","Comp_W_CO.W_CO1",
                                 "A3_Ind_D.1","Comp_W_CO.W_CO2")) %>%
  mutate(`Estimator Number`=row_number(),
         Type=c("GD","CPI","GD","GD","CPI","GD","GD","SA","SA","GD","SA"),
         Assumption=c("S4","S4","S2","S3","S3","S2","S2","S2","S2","S3","S2"),
         Name=c("GD_A4_j3","CPI_A4_j3","GD_A2_j3",
                "GD_A3_a2","CPI_A3_a2","GD_A2_a2",
                "GD_A2_a1","DCDH","CO1",
                "GD_A3_a1","CO2"))
TargetsSet_333 <- tibble(Estimator=c("A4_333_T.3","CPI.T_3","A2_333_T.3",
                                 "A3_333_D.2","CPI.D_2","A2_333_D.2",
                                 "A2_333_D.1","Comp_W_CH.W_M","Comp_W_CO.W_CO1",
                                 "A3_333_D.1","Comp_W_CO.W_CO2")) %>%
  mutate(`Estimator Number`=row_number(),
         Type=c("GD","CPI","GD","GD","CPI","GD","GD","SA","SA","GD","SA"),
         Assumption=c("S4","S4","S2","S3","S3","S2","S2","S2","S2","S3","S2"),
         Name=c("GD_A4_j3","CPI_A4_j3","GD_A2_j3",
                "GD_A3_a2","CPI_A3_a2","GD_A2_a2",
                "GD_A2_a1","DCDH","CO1",
                "GD_A3_a1","CO2"))

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

### Plots of Power:
Power13 <- 
  ggplot(Overall  %>% filter(SimNo %in% c(1:3)), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point() +
  facet_wrap(~SimNo, nrow=1, ncol=3) + theme_bw()

Power1 <- 
  ggplot(Overall  %>% filter(SimNo==1), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + theme_bw() +
  labs(x="Estimator",y="Type I Error (%)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(0,100), 
                     breaks=seq(0,100,by=20),
                     expand=c(0,2)) +
  geom_hline(yintercept=5, linetype="dashed", color="gray50")

Power2 <- 
  ggplot(Overall  %>% filter(SimNo==2), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + theme_bw() +
  labs(x="Estimator",y="Power (%)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(0,100), 
                     breaks=seq(0,100,by=20),
                     expand=c(0,2))

Power3 <- 
  ggplot(Overall  %>% filter(SimNo==3), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + theme_bw() +
  labs(x="Estimator",y="Power (%)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(0,100), 
                     breaks=seq(0,100,by=20),
                     expand=c(0,2))

Power46 <- 
  ggplot(Overall  %>% filter(SimNo %in% c(4:6), Assumption %in% c("S4","S5")), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point() +
  facet_wrap(~SimNo, nrow=1, ncol=3) + theme_bw()

Power4 <- 
  ggplot(Overall  %>% filter(SimNo==4, Assumption %in% c("S4","S5")), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + theme_bw() +
  labs(x="Estimator",y="Power (%)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(0,100), 
                     breaks=seq(0,100,by=20),
                     expand=c(0,2))

Power5 <- 
  ggplot(Overall  %>% filter(SimNo==5, Assumption %in% c("S4","S5")), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + theme_bw() +
  labs(x="Estimator",y="Power (%)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(0,100), 
                     breaks=seq(0,100,by=20),
                     expand=c(0,2))

Power6 <- 
  ggplot(Overall  %>% filter(SimNo==6, Assumption %in% c("S4","S5")), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + theme_bw() +
  labs(x="Estimator",y="Power (%)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(0,100), 
                     breaks=seq(0,100,by=20),
                     expand=c(0,2))

Power79 <- 
  ggplot(Overall  %>% filter(SimNo %in% c(7:9), Assumption %in% c("S3")), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point() +
  facet_wrap(~SimNo, nrow=1, ncol=3) + theme_bw()

Power7 <- 
  ggplot(Overall  %>% filter(SimNo==7, Assumption %in% c("S3")), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + theme_bw() +
  labs(x="Estimator",y="Power (%)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(0,100), 
                     breaks=seq(0,100,by=20),
                     expand=c(0,2))

Power8 <- 
  ggplot(Overall  %>% filter(SimNo==8, Assumption %in% c("S3")), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + theme_bw() +
  labs(x="Estimator",y="Power (%)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(0,100), 
                     breaks=seq(0,100,by=20),
                     expand=c(0,2))

Power9 <- 
  ggplot(Overall  %>% filter(SimNo==9, Assumption %in% c("S3")), 
         mapping=aes(x=`Estimator Number`, y=Power*100, 
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + theme_bw() +
  labs(x="Estimator",y="Power (%)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(0,100), 
                     breaks=seq(0,100,by=20),
                     expand=c(0,2))

### Plots of Estimates and SDs:
Est13 <- 
  ggplot(Overall %>% filter(SimNo %in% c(1:3)), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point() + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  facet_wrap(~SimNo, nrow=1, ncol=3)

Est1 <- 
  ggplot(Overall %>% filter(SimNo==1), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray50")

Est2 <- 
  ggplot(Overall %>% filter(SimNo==2), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_hline(yintercept=-0.02, linetype="dashed", color="gray50")

Est3 <- 
  ggplot(Overall %>% filter(SimNo==3), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point(size=2) + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_hline(yintercept=-0.04, linetype="dashed", color="gray50")

Est46 <- 
  ggplot(Overall %>% filter(SimNo %in% c(4:6),
                            Assumption %in% c("S2","S4","S5")), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point() + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  facet_wrap(~SimNo, nrow=1, ncol=3)

Est4 <- 
  ggplot(Overall %>% filter(SimNo==4,
                            Assumption %in% c("S2","S4","S5")), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point() + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.02, x=0.5, xend=7.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=0.005, x=11.5, xend=13.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.0033333, x=13.5, xend=16.5,
               linetype="dashed", color="grey50")

Est5 <- 
  ggplot(Overall %>% filter(SimNo==5,
                            Assumption %in% c("S2","S4","S5")), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point() + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.02, x=0.5, xend=7.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=0.005694444, x=11.5, xend=13.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.001904762, x=13.5, xend=16.5,
               linetype="dashed", color="grey50")
             
Est6 <- 
  ggplot(Overall %>% filter(SimNo==6,
                            Assumption %in% c("S2","S4","S5")), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point() + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.02, x=0.5, xend=7.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.0105, x=11.5, xend=13.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.01428571, x=13.5, xend=16.5,
               linetype="dashed", color="grey50")

Est79 <- 
  ggplot(Overall %>% filter(SimNo %in% c(7:9),
                            Assumption %in% c("S2","S3","S5")), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point() + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  facet_wrap(~SimNo, nrow=1, ncol=3)

Est7 <- 
  ggplot(Overall %>% filter(SimNo==7,
                            Assumption %in% c("S2","S3","S5")), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point() + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.025, x=7.5, xend=9.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.0225, x=9.5, xend=11.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-.01625, x=11.5, xend=13.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.01833333, x=13.5, xend=16.5,
               linetype="dashed", color="grey50")

Est8 <- 
  ggplot(Overall %>% filter(SimNo==8,
                            Assumption %in% c("S2","S3","S5")), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point() + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.0214, x=7.5, xend=9.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.02, x=9.5, xend=11.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.0105, x=11.5, xend=13.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.01428571, x=13.5, xend=16.5,
               linetype="dashed", color="grey50")

Est9 <- 
  ggplot(Overall %>% filter(SimNo==9,
                            Assumption %in% c("S2","S3","S5")), 
         mapping=aes(x=`Estimator Number`, y=`Mean Estimate`, 
                     ymin=Lower, ymax=Upper,
                     color=Type, shape=Assumption)) + 
  geom_point() + geom_errorbar() +
  theme_bw() + labs(x="Estimator", y="Estimate (Mean \U00B1 SD)") +
  scale_x_discrete() + 
  scale_y_continuous(limits=c(-0.10,0.015), 
                     breaks=seq(-0.1,0.01,by=0.02),
                     expand=c(0,0)) +
  geom_segment(y=-0.01, x=7.5, xend=9.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.02, x=9.5, xend=11.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.045, x=11.5, xend=13.5,
               linetype="dashed", color="grey50") +
  geom_segment(y=-0.036666667, x=13.5, xend=16.5,
               linetype="dashed", color="grey50")

### Export Plots
ggsave(filename=paste0(outdir,"Sim_Power_1.png"),
       plot=Power1 + guides(color="none", shape="none") + labs(x=NULL, title="A) Scenario 1") +
         Power2 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="B) Scenario 2") +
         Power3 + guides(color="none") +labs(x=NULL, y=NULL, title="C) Scenario 3") +
         Power4 + guides(color="none", shape="none") + labs(x=NULL, title="D) Scenario 4") +
         Power5 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="E) Scenario 5") +
         Power6 + guides(shape="none") + labs(x=NULL, y=NULL, title="F) Scenario 6") +
         Power7 + guides(color="none", shape="none") + labs(title="G) Scenario 7") +
         Power8 + guides(color="none", shape="none") + labs(y=NULL, title="H) Scenario 8") +
         Power9 + guides(color="none", shape="none") + labs(y=NULL, title="I) Scenario 9") +
         plot_layout(nrow=3, ncol=3, byrow=TRUE),
       width=8, height=7, units="in", dpi=300)
ggsave(filename=paste0(outdir,"Sim_Ests_1.png"),
       plot=Est1 + guides(color="none", shape="none") + labs(x=NULL, title="A) Scenario 1") +
         Est2 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="B) Scenario 2") +
         Est3 + guides(color="none") + labs(x=NULL, y=NULL, title="C) Scenario 3") +
         Est4 + guides(color="none", shape="none") + labs(x=NULL, title="D) Scenario 4") +
         Est5 + guides(color="none", shape="none") + labs(x=NULL, y=NULL, title="E) Scenario 5") +
         Est6 + guides(shape="none") + labs(x=NULL, y=NULL, title="F) Scenario 6") +
         Est7 + guides(color="none", shape="none") + labs(title="G) Scenario 7") +
         Est8 + guides(color="none", shape="none") + labs(y=NULL, title="H) Scenario 8") +
         Est9 + guides(color="none", shape="none") + labs(y=NULL, title="I) Scenario 9") +
         plot_layout(nrow=3, ncol=3, byrow=TRUE),
       width=8, height=7, units="in", dpi=300)
