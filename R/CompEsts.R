#######################################
###### File: CompEsts.R ###############
###### Lee Kennedy-Shaffer ############
###### Created 2024/07/18 #############
#######################################

## Load packages
require(dplyr)
require(tibble)
require(tidyr)
require(bbmle)

## From a DFT output object,
### derive the weights for various other estimators
Comp_Ests <- function(DFT_obj,
                      estimator=c("TW","CS","SA","CH","CO")) {
  ### First, augment D with start times:
  Starts <- DFT_obj$Theta$Full %>%
    dplyr::select(Cl.Num,Start_j) %>%
    distinct()
  if (max(Starts$Cl.Num) < DFT_obj$N) {
    Starts <- Starts %>%
      bind_rows(tibble(Cl.Num=(max(Starts$Cl.Num)+1):DFT_obj$N,
                                   Start_j=Inf))
  }
  D_use <- DFT_obj$D_aug %>% dplyr::select(i,i.prime,j,j.prime,Type,TypeLabel) %>%
    left_join(Starts %>% rename(i.start=Start_j),
              by=join_by(i == Cl.Num)) %>%
    left_join(Starts %>% rename(i.prime.start=Start_j),
              by=join_by(i.prime == Cl.Num))

  res <- NULL

  if ("TW" %in% estimator) { ## Based on weights from Goodman-Bacon (2021):
    Starts_n <- Starts %>% group_by(Start_j) %>%
      dplyr::summarize(n.Clust=n()) %>%
      dplyr::mutate(n=n.Clust/sum(n.Clust),
                    Dbar=if_else(is.infinite(Start_j),0,
                            (DFT_obj$J-Start_j+1)/DFT_obj$J)) %>%
      dplyr::rename(G=Start_j)
    Cross_starts <- cross_join(Starts_n, Starts_n,
                               suffix=c(".k",".l")) %>%
      dplyr::filter(G.k != G.l, Dbar.k > 0)
    TW.Weights <- Cross_starts %>%
      dplyr::mutate(s.num=if_else(Dbar.l==0, n.k*(n.l)*Dbar.k*(1-Dbar.k),
                                  if_else(G.k < G.l, n.k*n.l*(Dbar.k-Dbar.l)*(1-Dbar.k),
                                          -1*n.k*n.l*Dbar.k*(Dbar.l-Dbar.k))),
                    n.DID=if_else(Dbar.l==0, (DFT_obj$J-G.k+1)*(G.k-1)*n.Clust.k*n.Clust.l,
                                  if_else(G.k < G.l, (G.l-G.k)*(G.k-1)*n.Clust.k*n.Clust.l,
                                          (DFT_obj$J-G.k+1)*(G.k-G.l)*n.Clust.k*n.Clust.l))) %>%
      dplyr::mutate(s=s.num/sum(abs(s.num)),
                    W_TW=s/n.DID)
    TW_w <- D_use %>% mutate(G.k=if_else(Type==2, i.start,
                                         if_else(Type==5, i.prime.start, NA)),
                              G.l=if_else(Type==2, i.prime.start,
                                          if_else(Type==5, i.start, NA))) %>%
      left_join(TW.Weights %>% dplyr::select(G.k,G.l,W_TW),
                by=c("G.k","G.l")) %>%
      dplyr::mutate(W_TW=if_else(is.na(W_TW), 0, W_TW))
    res <- c(res, list(W_TW=TW_w %>% dplyr::select(W_TW)))
  }

  if ("CS" %in% estimator) {
    Group <- unique(Starts$Start_j[Starts$Start_j < max(Starts$Start_j)])
    Time <- unique(D_use %>% filter(Type %in% 1:3) %>% pull(j.prime))
    CS_key <- crossing(Group,Time) %>% mutate(Col=1:n())
    CS_gt_w_int <- apply(as.matrix(CS_key),
                          MARGIN=1,
                          FUN=function(row) CS_wt_fun(D_use, row["Group"], row["Time"]))
    CS_gt_w <- apply(CS_gt_w_int, 2, FUN=function(x) x/sum(abs(x)))
    CS_key_agg <- CS_key %>% mutate(Type=if_else(Time >= Group,"Post","Pre"),
                                      EventTime=Time-Group) %>%
      left_join(Starts %>% group_by(Start_j) %>% dplyr::summarize(Num=n()),
                by=join_by(Group==Start_j))
    CS_key_agg <- CS_key_agg %>%
      left_join(CS_key_agg %>% group_by(EventTime) %>% dplyr::summarize(Num_dynamic=sum(Num))) %>%
      left_join(CS_key_agg %>% group_by(Group) %>%
                  dplyr::summarize(Num_group=sum(if_else(Type=="Post",Num,0)))) %>%
      left_join(CS_key_agg %>% group_by(Time) %>%
                  dplyr::summarize(Num_calendar=sum(if_else(Type=="Post",Num,0)))) %>%
      mutate(W_simple=if_else(Type=="Post",Num,0)/sum(if_else(Type=="Post",Num,0)),
             W_dynamic=(if_else(Type=="Post",Num/Num_dynamic,0)/sum(if_else(Type=="Post",Num/Num_dynamic,0))),
             W_group=(if_else(Type=="Post",Num/Num_group,0)/sum(if_else(Type=="Post",Num/Num_group,0))),
             W_calendar=(if_else(Type=="Post",Num/Num_calendar,0)/sum(if_else(Type=="Post",Num/Num_calendar,0))))
    CS_w <- CS_gt_w %*% as.matrix(CS_key_agg %>% dplyr::select(starts_with("W_")))
    res <- c(res, list(Key_CS=CS_key_agg, W_CS_gt=CS_gt_w, W_CS=CS_w))
  }

  if ("SA" %in% estimator) {
    Rel_Pds <- unique(D_use %>% filter(Type %in% 1:3) %>%
                        mutate(RelPd=j.prime-i.start) %>% pull(RelPd))
    Rel_Pds <- sort(Rel_Pds[Rel_Pds != -1])
    Ctrl.is <- Starts %>% filter(Start_j==max(Start_j)) %>% pull(Cl.Num)
    CATT_el_key <- crossing(Period=Rel_Pds,
                            Cohort=Starts %>% filter(Start_j != max(Start_j)) %>%
                              pull(Start_j))
    CATT_el_key$Num <- apply(as.matrix(CATT_el_key), MARGIN=1,
                             FUN=function(row) nrow(D_use %>% filter(Type %in% 1:3,
                                                                    i.start==row["Cohort"],
                                                                    i.prime %in% Ctrl.is,
                                                                    j.prime-i.start==row["Period"],
                                                                    j-i.start==-1)) +
                               nrow(D_use %>% filter(Type %in% 1:3,
                                                    i.start==row["Cohort"],
                                                    i.prime %in% Ctrl.is,
                                                    j.prime-i.start==-1,
                                                    j-i.start==row["Period"])))
    CATT_el_key <- CATT_el_key %>% filter(Num != 0) %>%
      mutate(Column=1:n())
    CATT_el_w_int <- apply(CATT_el_key, MARGIN=1,
                           FUN=function(row) SA_wt_fun(D_use, Ctrl.is,
                                                       row["Period"],row["Cohort"]))
    CATT_el_w <- apply(CATT_el_w_int, MARGIN=2,
                       FUN=function(x) x/sum(abs(x)))
    SA_l_key <- CATT_el_key %>% dplyr::select(Period,Num) %>% group_by(Period) %>%
      dplyr::summarize(Num=sum(Num)) %>% mutate(Column=1:n()) %>%
      mutate(Type=if_else(Period >= 0,"Post","Pre"),
             W_SA.W_ATT=if_else(Type=="Post",Num,0)/sum(if_else(Type=="Post",Num,0)))
    SA_l_w <- apply(SA_l_key,
                    MARGIN=1,
                    FUN=function(row) SA_CATT_agg(CATT_el_key,
                                               CATT_el_w,
                                               as.numeric(row["Period"]))/as.numeric(row["Num"]))
    res <- c(res,
             list(Key_SA_l=SA_l_key, W_SA_l=SA_l_w,
             W_SA=SA_l_w %*% as.matrix(SA_l_key %>% dplyr::select(W_SA.W_ATT))))
  }

  if ("CH" %in% estimator) {
    CH_Pt_key <- D_use %>% dplyr::filter(Type==2,
                                         j.prime-j==1,
                                         j.prime==i.start) %>%
      group_by(j.prime) %>%
      dplyr::summarize(Num=n(),
                       Num.Trt=length(unique(i))) %>%
      rename(Period=j.prime) %>%
      mutate(Column=1:n())
    CH_Pt_w_int <- apply(CH_Pt_key, MARGIN=1,
                           FUN=function(row) CH_wt_fun(D_use,row["Period"]))
    CH_Pt_w <- apply(CH_Pt_w_int, 2, function(x) x/sum(abs(x)))
    res <- c(res, list(W_CH=tibble(W_CH.W_M = CH_Pt_w %*% CH_Pt_key$Num.Trt/sum(CH_Pt_key$Num.Trt))))
  }

  if ("CO" %in% estimator) {
    CO_Pt_key <- D_use %>% dplyr::filter(Type %in% c(2,5),
                                         j.prime-j==1,
                                         j.prime==i.start | j.prime==i.prime.start) %>%
      group_by(j.prime,Type) %>%
      dplyr::summarize(Num.Comps=n(),
                       Num.i=length(unique(i)),
                       Num.i.prime=length(unique(i.prime))) %>%
      pivot_wider(id_cols="j.prime",
                  names_from="Type",
                  values_from=starts_with("Num"),
                  names_sep=".",
                  values_fill=0) %>%
      mutate(Num.Comps=Num.Comps.2+Num.Comps.5,
             Num.Trt=max(Num.i.2,Num.i.prime.5),
             Num.Ctrl=Num.i.prime.2+Num.i.5,
             Num.Comps.Res=Num.Comps.2,
             Num.Trt.Res=Num.i.2,
             Num.Ctrl.Res=Num.i.prime.2) %>%
      ungroup() %>%
      dplyr::select(j.prime,Num.Comps,Num.Trt,Num.Ctrl,
                    Num.Comps.Res,Num.Trt.Res,Num.Ctrl.Res) %>%
      rename(Period=j.prime) %>%
      mutate(Column=1:n())
    CO_Pt_w_int <- apply(CO_Pt_key, MARGIN=1,
                            FUN=function(row) CO_wt_fun(D_use,row["Period"]))
    CO_Pt_All <- apply(CO_Pt_w_int, 2,
                      function(x) x/if_else(sum(abs(x))==0,1,sum(abs(x))))
    CO_Pt_All <- apply(CO_Pt_w_int,
                       2, function(x) x/sum(abs(x)))
    CO_Pt_Res <- apply(CO_Pt_w_int,
                       2, function(x) if_else(x<0,0,x)/if_else(sum(x[x>=0])==0,1,sum(x[x>=0])))
    CO_Pt_key <- CO_Pt_key %>%
      mutate(W_CO1=if_else(Num.Comps.Res > 0,1,0)/sum(if_else(Num.Comps.Res > 0,1,0)),
             W_CO2=if_else(Num.Comps.Res > 0,
                           2/(1/Num.Trt.Res+1/Num.Ctrl.Res),
                           0)/sum(if_else(Num.Comps.Res > 0,
                                          2/(1/Num.Trt.Res+1/Num.Ctrl.Res),
                                          0)),
             W_CO3=if_else(Num.Comps > 0,1,0)/sum(if_else(Num.Comps > 0,1,0)))
    CO_w <- bind_cols(CO_Pt_Res %*% as.matrix(CO_Pt_key %>% dplyr::select(c("W_CO1","W_CO2"))),
                      CO_Pt_All %*% as.matrix(CO_Pt_key %>% dplyr::select("W_CO3")))
    res <- c(res,list(Key_CO=CO_Pt_key,
                      W_CO_Pt_All=CO_Pt_All,
                      W_CO_Pt_Res=CO_Pt_Res,
                      W_CO=CO_w))
  }
 Weights <- do.call(cbind, res[paste0("W_",estimator)])
  return(list(Weights=Weights,Full=res))
}

Comp_Ests_Weights <- function(DFT_obj, Amat,
                              estimator=c("TW","CS","SA","CH","CO","NP")) {
  DID.weights <- NULL
  Obs.weights <- NULL
  if (sum(c("TW","CS","SA","CH","CO") %in% estimator) > 0) {
    DID.ests <- estimator[estimator %in% c("TW","CS","SA","CH","CO")]
    DID.weights <- bind_cols(DID.weights,Comp_Ests(DFT_obj,
                                                   DID.ests)$Weights)
    Obs.weights <- bind_cols(Obs.weights, t(Amat) %*% as.matrix(DID.weights))
  }
  if ("NP" %in% estimator) {
    N_Trt <- apply(DFT_obj$Theta$Schematic, 2, function(x) sum(x > 0))
    N_Ctrl <- DFT_obj$N - N_Trt
    W_Trt <- if_else(N_Trt==0 | N_Ctrl==0,0,1/N_Trt)
    W_Ctrl <- if_else(N_Trt==0 | N_Ctrl==0,0,-1/(N_Ctrl))
    Mat <- matrix(c(W_Trt, N_Trt, W_Ctrl, N_Ctrl), nrow=4, byrow=TRUE)
    W_Mat <- apply(Mat, 2, FUN=function(col) c(rep(col[1], col[2]), rep(col[3], col[4])))
    NP_w_int <- as.vector(t(W_Mat))
    Mult_Eq <- if_else(N_Trt > 0 & N_Ctrl > 0, 1, 0)/sum(if_else(N_Trt > 0 & N_Ctrl > 0, 1, 0))
    Mult_ATT <- if_else(N_Trt > 0 & N_Ctrl > 0, N_Trt, 0)/sum(if_else(N_Trt > 0 & N_Ctrl > 0, N_Trt, 0))
    Mult_IV <- if_else(N_Trt > 0 & N_Ctrl > 0, (1/N_Trt + 1/N_Ctrl)^(-1),0)/sum(if_else(N_Trt > 0 & N_Ctrl > 0, (1/N_Trt + 1/N_Ctrl)^(-1),0))
    Obs.weights <- bind_cols(Obs.weights,
                             W_NP_Eq=NP_w_int * rep(Mult_Eq, times=DFT_obj$N),
                             W_NP_ATT=NP_w_int * rep(Mult_ATT, times=DFT_obj$N),
                             W_NP_IV=NP_w_int * rep(Mult_IV, times=DFT_obj$N))
  }

  return(list(DID.weights=DID.weights,
              Obs.weights=Obs.weights))
}

### Helper functions:
CS_wt_fun <- function(D_use,group,time) {
  if (time >= group) {
    if_else(D_use$Type %in% 1:3 & D_use$i.start==group & D_use$i.prime.start != group &
              D_use$j.prime==time & D_use$j==group-1,
            1,0)
    # apply(D_use, MARGIN=1,
    #       FUN=function(r) if_else(r["Type"] %in% 1:3 & r["i.start"]==group &
    #                                 r["i.prime.start"] != group & r["j.prime"]==time &
    #                                 r["j"]==group-1,1,0))
  } else {
    if_else(D_use$Type & D_use$j.prime & D_use$j==time-1,
            if_else(D_use$i.start==group & D_use$i.prime.start != group, 1,
                    if_else(D_use$i.start != group & D_use$i.prime.start == group, -1, 0)), 0)
    # apply(D_use, MARGIN=1,
    #       FUN=function(r) if_else(r["Type"]==1 & r["j.prime"]==time & r["j"]==time-1,
    #                               if_else(r["i.start"]==group & r["i.prime.start"] != group,1,
    #                                       if_else(r["i.start"] != group & r["i.prime.start"] == group,-1,0)),0))
  }
}

SA_wt_fun <- function(D_use,Ctrl.is,period,cohort) {
  if (period >= 0) {
    if_else(D_use$Type %in% 1:3 & D_use$i.start==cohort &
              D_use$i.prime %in% Ctrl.is &
              as.numeric(D_use$j.prime) - as.numeric(D_use$i.start) == period &
              as.numeric(D_use$j) - as.numeric(D_use$i.start) == -1,
            1, 0)
    # apply(D_use, MARGIN=1,
    #       FUN=function(r) if_else(r["Type"] %in% 1:3 &
    #                                 r["i.start"]==cohort &
    #                                 r["i.prime"] %in% Ctrl.is &
    #                                 as.numeric(r["j.prime"])-as.numeric(r["i.start"])==period &
    #                                 as.numeric(r["j"])-as.numeric(r["i.start"])==-1,
    #                               1,0))
  } else {
    if_else(D_use$Type %in% 1:3 & D_use$i.start==cohort &
              D_use$i.prime %in% Ctrl.is &
              as.numeric(D_use$j.prime) - as.numeric(D_use$i.start) == -1 &
              as.numeric(D_use$j) - as.numeric(D_use$i.start) == period,
            -1, 0)
    # apply(D_use, MARGIN=1,
    #       FUN=function(r) if_else(r["Type"] %in% 1:3 &
    #                                 r["i.start"]==cohort &
    #                                 r["i.prime"] %in% Ctrl.is &
    #                                 as.numeric(r["j.prime"])-as.numeric(r["i.start"])==-1 &
    #                                 as.numeric(r["j"])-as.numeric(r["i.start"])==period,
    #                               -1,0))
  }
}

SA_CATT_agg <- function(key_tbl, wt_tbl, period) {
  key_vals <- key_tbl %>% dplyr::filter(Period==period)
  wt_tbl[,key_vals$Column,drop=FALSE] %*% key_vals$Num
}

CH_wt_fun <- function(D_use,period) {
  if_else(D_use$Type==2 & D_use$j.prime==D_use$i.start &
            D_use$j.prime==period & D_use$j == period-1,
          1,0)
  # apply(D_use, MARGIN=1,
  #       FUN=function(r) if_else(r["Type"]==2 &
  #                                 r["j.prime"]==r["i.start"] &
  #                                 r["j.prime"]==period &
  #                                 r["j"]==period-1,
  #                               1,0))
}

CO_wt_fun <- function(D_use,period) {
  if_else(D_use$Type==2,
          if_else(D_use$j.prime==D_use$i.start &
                    D_use$j.prime==period &
                    D_use$j==period-1, 1, 0),
          if_else(D_use$Type==5,
                  if_else(D_use$j.prime == D_use$i.prime.start &
                            D_use$j.prime == period & D_use$j == period-1,
                          -1,0),0))
  # apply(D_use, MARGIN=1,
  #       FUN=function(r) if_else(r["Type"]==2,
  #                               if_else(r["j.prime"]==r["i.start"] &
  #                                         r["j.prime"]==period &
  #                                         r["j"]==period-1,1,0),
  #                               if_else(r["Type"]==5,
  #                                       if_else(r["j.prime"]==r["i.prime.start"] &
  #                                                 r["j.prime"]==period &
  #                                                 r["j"]==period-1,-1,0),0)))
}


#### CLWP and CLWPA analysis functions ####
## See Voldal et al., Statist Med, 2024, https://doi.org/10.1002/sim.10120 ##
## and its supplemental material for details on these methods and the ##
## underlying code. If using these methods, please refer to and cite ##
## Voldal et al. 2024 ##

#####
#Helper functions for the unadjusted vertical CL
#####

#Composite likelihood
#Inputs:
#theta - a value for the treatment effect
#sigma - a value for the standard deviation of the vertical contrasts
#my.data.full - a dataset in the format produced by swCRTdesign
#Returns: negative log likelihood
loglik_vertical_CML <- function(theta, sigma, my.data.full){
  #First, we need to go from a data frame of cluster-period observations to a data frame with
  #every vertical contrast between treatment and control cluster-periods.
  my.data.combos <- data.frame(treated=double(),control=double())
  for(j in unique(my.data.full$Period)){
    treated_list.j <- my.data.full$Y.ij.bar[my.data.full$Period == j & my.data.full$Interv == 1]
    control_list.j <- my.data.full$Y.ij.bar[my.data.full$Period == j & my.data.full$Interv == 0]
    combos.j <- expand.grid(treated=treated_list.j,control=control_list.j)
    my.data.combos <- rbind(my.data.combos,combos.j)
  }
  my.data.combos$difference <- my.data.combos$treated-my.data.combos$control
  R=suppressWarnings(dnorm(x=my.data.combos$difference,mean=theta,sd=sigma))
  return(-sum(log(R)))
}

#Inputs:
#fitted_params - fitted parameters from the vertical CL, e.g. my.fitted.cl@coef, a vector of (theta, sigma)
#difference - one vertical difference
#Returns: u (composite score function) for one vertical difference
get_u_per_contrast <- function(fitted_params,difference){
  theta_cl <- fitted_params["theta"]
  sigmasq_cl <- fitted_params["sigma"]^2#Note: going from sd to var

  u_1 <- (1/sigmasq_cl)*(difference-theta_cl)

  u_2 <- -(1/2)*(1/sigmasq_cl)+(1/2)*(1/sigmasq_cl^2)*(difference-theta_cl)^2

  return(c(u_1,u_2))

}

#Inputs:
#fitted_params - fitted parameters from the vertical CL, e.g. my.fitted.cl@coef, a vector of (theta, sigma)
#difference - one vertical difference
#Returns: a matrix of the second derivatives (derivative of the score) for one vertical difference
get_triangledown_u_per_contrast <- function(fitted_params,difference){
  theta_cl <- fitted_params["theta"]
  sigmasq_cl <- fitted_params["sigma"]^2#Note: going from sd to var

  block_1_1 <- -1/sigmasq_cl

  block_1_2 <- -1/sigmasq_cl^2*(difference-theta_cl)

  block_2_2 <- (1/2)*(1/sigmasq_cl^2)-(1/sigmasq_cl^3)*(difference-theta_cl)^2

  row_1 <- c(block_1_1,block_1_2)
  row_2 <- c(block_1_2,block_2_2)

  full_matrix <- rbind(row_1,row_2)

  return(full_matrix)

}


## ---------------------------------------------------------------------------------
#####
#Helper functions for the adjusted vertical CL
#####

#Composite likelihood
#Inputs:
#theta - a value for the treatment effect
#sigma - a value for the standard deviation of the vertical contrasts, after adjusting for baseline
#betabase - a value for the coefficient for the baseline difference
#my.data.full - a dataset in the format produced by swCRTdesign
#Returns: negative log likelihood
loglik_vertical_CML_adj <- function(theta, sigma, betabase, my.data.full){
  #First, we need to go from a data frame of cluster-period observations to a data frame with
  #every vertical contrast between treatment and control cluster-periods.
  my.data.combos <- data.frame(treated=double(),control=double())
  my.data.base.combos <- data.frame(treated=double(),control=double())
  for(j in unique(my.data.full$Period)){
    treated_list.j <- my.data.full$Y.ij.bar[my.data.full$Period == j & my.data.full$Interv == 1]
    control_list.j <- my.data.full$Y.ij.bar[my.data.full$Period == j & my.data.full$Interv == 0]
    combos.j <- expand.grid(treated=treated_list.j,control=control_list.j)
    my.data.combos <- rbind(my.data.combos,combos.j)
    #Now do the same thing in the same order, but get the baseline measurements
    #Note: assuming here that the first period is all-control
    #(so won't be contributing to vertical contrasts, and don't need to account for differing trt status)
    treated_list_base.j <- my.data.full$Y.ij.bar[my.data.full$Period == 1 &
                                                       my.data.full$Cluster %in% my.data.full$Cluster[my.data.full$Period == j & my.data.full$Interv == 1]]
    control_list_base.j <- my.data.full$Y.ij.bar[my.data.full$Period == 1 &
                                                       my.data.full$Cluster %in% my.data.full$Cluster[my.data.full$Period == j & my.data.full$Interv == 0]]
    combos.base.j <- expand.grid(treated=treated_list_base.j,control=control_list_base.j)
    my.data.base.combos <- rbind(my.data.base.combos,combos.base.j)
  }
  my.data.combos$difference <- my.data.combos$treated-my.data.combos$control
  my.data.base.combos$difference <- my.data.base.combos$treated-my.data.base.combos$control
  R=suppressWarnings(dnorm(x=my.data.combos$difference,mean=theta+betabase*my.data.base.combos$difference,sd=sigma))
  return(-sum(log(R)))
}


#Inputs:
#fitted_params - fitted parameters from the vertical CL, e.g. my.fitted.cl@coef, a vector of (theta, sigma, betabase)
#difference - one vertical difference
#difference_base - the analogous difference at baseline
#Returns: u (composite score function) for one vertical difference
get_u_per_contrast_adj <- function(fitted_params,difference,difference_base){
  theta_cl <- fitted_params["theta"]
  sigmasq_cl <- fitted_params["sigma"]^2#Note: going from sd to var
  betabase_cl <- fitted_params["betabase"]

  u_1 <- (1/sigmasq_cl)*(difference-(theta_cl+betabase_cl*difference_base))

  u_2 <- -(1/2)*(1/sigmasq_cl)+(1/2)*(1/sigmasq_cl^2)*(difference-(theta_cl+betabase_cl*difference_base))^2

  u_3 <- (1/sigmasq_cl)*(difference-(theta_cl+betabase_cl*difference_base))*difference_base


  return(c(u_1,u_2,u_3))

}

#Inputs:
#fitted_params - fitted parameters from the vertical CL, e.g. my.fitted.cl@coef, a vector of (theta, sigma, betabase)
#difference - one vertical difference
#difference_base - the analogous difference at baseline
#Returns: a matrix of the second derivatives (derivative of the score) for one vertical difference
get_triangledown_u_per_contrast_adj <- function(fitted_params,difference,difference_base){
  theta_cl <- fitted_params["theta"]
  sigmasq_cl <- fitted_params["sigma"]^2#Note: going from sd to var
  betabase_cl <- fitted_params["betabase"]

  block_1_1 <- -1/sigmasq_cl

  block_1_2 <- -1/sigmasq_cl^2*(difference-(theta_cl+betabase_cl*difference_base))

  block_1_3 <- -difference_base/sigmasq_cl


  block_2_2 <- (1/2)*(1/sigmasq_cl^2)-(1/sigmasq_cl^3)*(difference-(theta_cl+betabase_cl*difference_base))^2

  block_2_3 <- -difference_base/sigmasq_cl^2*(difference-(theta_cl+betabase_cl*difference_base))

  block_3_3 <- -difference_base^2/sigmasq_cl

  row_1 <- c(block_1_1,block_1_2,block_1_3)
  row_2 <- c(block_1_2,block_2_2,block_2_3)
  row_3 <- c(block_1_3,block_2_3,block_3_3)

  full_matrix <- rbind(row_1,row_2,row_3)

  return(full_matrix)

}

#####
#Composite likelihood within period (CLWP)

#using starting values from GEE
#Note - I did do some experimenting, and it doesn't seem to be very sensitive to starting values.
# start_theta <- summary(gee.independence)$coef["tx.var","Estimate"]
# start_sigma <- sqrt(2*var(gee.independence$residuals))
CLWP_fit <- function(my.data, start_theta, start_sigma, N) {
  cml.vertical <-mle2(function(theta,sigma){loglik_vertical_CML(theta=theta, sigma=sigma,my.data.full=my.data)},
                      start=list(theta=start_theta, sigma=start_sigma))

  #Recreate data frame of all the differences
  my.data.combos <- data.frame(treated=double(),control=double(),cluster_trt=double(),cluster_ctrl=double())
  for(j in unique(my.data$Period)){
    treated_list.j <- my.data$Y.ij.bar[my.data$Period == j & my.data$Interv == 1]
    control_list.j <- my.data$Y.ij.bar[my.data$Period == j & my.data$Interv == 0]
    combos.j <- expand.grid(treated=treated_list.j,control=control_list.j)
    #get cluster ID's for each contrast to pull later
    cluster_treated_list.j <- my.data$Cluster[my.data$Period == j & my.data$Interv == 1]
    cluster_control_list.j <- my.data$Cluster[my.data$Period == j & my.data$Interv == 0]
    cluster_combos.j <- expand.grid(cluster_trt=cluster_treated_list.j,cluster_ctrl=cluster_control_list.j)
    my.data.combos <- rbind(my.data.combos,cbind(combos.j,cluster_combos.j))
  }
  my.data.combos$difference <- my.data.combos$treated-my.data.combos$control


  #Set up empty matrices to sum
  triangledown_u_sum <- matrix(0,nrow=2,ncol=2)#2x2 for the two parameters
  u_sq_sum <- matrix(0,nrow=2,ncol=2)


  for(i in 1:nrow(my.data.combos)){
    triangledown_u_sum <- triangledown_u_sum +get_triangledown_u_per_contrast(fitted_params=cml.vertical@coef,
                                                                              difference=my.data.combos$difference[i])
  }

  #To do J, need to sum u's within clusters, then square, THEN sum over clusters
  for(i in unique(my.data$Cluster)){
    u.i <- c(0,0)#initializing score vector (1x2)
    contrast_involves_i <- my.data.combos$cluster_trt == i | my.data.combos$cluster_ctrl == i
    my.data.combos.i <- my.data.combos[contrast_involves_i,]
    for(k in 1:nrow(my.data.combos.i)){#loop through all the contrasts that involve this cluster
      u.ik <- get_u_per_contrast(fitted_params=cml.vertical@coef,difference=my.data.combos.i$difference[k])
      u.i <- u.i+u.ik
    }
    u_sq_sum <- u_sq_sum + u.i %*% t(u.i)
  }


  my_H <- - triangledown_u_sum/N
  my_J <- u_sq_sum/N
  my_G <- my_H %*% solve(my_J) %*% my_H
  my_variance_matrix <- solve(my_G)

  return(c(CLWP_est=cml.vertical@coef["theta"],
           CLWP_se=((1/N)*my_variance_matrix[1,1])^.5,
           CLWP_pval=2 * pnorm(abs((cml.vertical@coef["theta"])/(((1/N)*my_variance_matrix[1,1])^.5)),
                               lower.tail = FALSE)))
}


CLWPA_fit <- function(my.data, start_theta, start_sigma, N) {

  #using same starting values from GEE as the unadjusted version
  #the starting standard deviation will be a little too large for the adjusted version
  #based on informal experimenting, starting vals too big are fine (too small are riskier).
  #If desired, could fit a GEE that adjusted for baseline to get a better start.
  cml.vertical.adj <- mle2(function(theta,sigma,betabase){loglik_vertical_CML_adj(theta=theta, sigma=sigma,betabase=betabase,
                                                                                  my.data.full=my.data)},start=list(theta=start_theta, sigma=start_sigma,betabase=0))
  #Recreate data frame of all the differences
  #####
  my.data.combos <- data.frame(treated=double(),control=double(),cluster_trt=double(),cluster_ctrl=double())
  for(j in unique(my.data$Period)){
    treated_list.j <- my.data$Y.ij.bar[my.data$Period == j & my.data$Interv == 1]
    control_list.j <- my.data$Y.ij.bar[my.data$Period == j & my.data$Interv == 0]
    combos.j <- expand.grid(treated=treated_list.j,control=control_list.j)
    cluster_treated_list.j <- my.data$Cluster[my.data$Period == j & my.data$Interv == 1]
    cluster_control_list.j <- my.data$Cluster[my.data$Period == j & my.data$Interv == 0]
    cluster_combos.j <- expand.grid(cluster_trt=cluster_treated_list.j,cluster_ctrl=cluster_control_list.j)
    my.data.combos <- rbind(my.data.combos,cbind(combos.j,cluster_combos.j))
  }
  my.data.combos$difference <- my.data.combos$treated-my.data.combos$control

  #use cluster ID's in cluster_trt and cluster_ctrl to fill in baseline response values
  #note: merge re-organizes order; that's fine for the CL function
  my.data.combos <- merge(my.data.combos,my.data[my.data$Period == 1,
                                                 c("Cluster","Y.ij.bar")],by.x="cluster_trt",by.y="Cluster",all.x=TRUE)
  my.data.combos$treated_base <- my.data.combos$Y.ij.bar
  my.data.combos <- merge(my.data.combos,my.data[my.data$Period == 1,
                                                 c("Cluster","Y.ij.bar")],by.x="cluster_ctrl",by.y="Cluster",all.x=TRUE)
  my.data.combos$control_base <- my.data.combos$Y.ij.bar.y

  #difference in baseline
  my.data.combos$difference_base <- my.data.combos$treated_base-my.data.combos$control_base


  #Set up empty matrices to sum
  triangledown_u_sum <- matrix(0,nrow=3,ncol=3)
  u_sq_sum <- matrix(0,nrow=3,ncol=3)



  for(i in 1:nrow(my.data.combos)){
    triangledown_u_sum <- triangledown_u_sum +get_triangledown_u_per_contrast_adj(fitted_params=cml.vertical.adj@coef,
                                                                                  difference=my.data.combos$difference[i],difference_base=my.data.combos$difference_base[i])
  }

  for(i in unique(my.data$Cluster)){
    u.i <- c(0,0,0)#initializing score vector (1x3)
    contrast_involves_i <- my.data.combos$cluster_trt == i | my.data.combos$cluster_ctrl == i
    my.data.combos.i <- my.data.combos[contrast_involves_i,]
    for(k in 1:nrow(my.data.combos.i)){#loop through all the contrasts that involve this cluster
      u.ik <- get_u_per_contrast_adj(fitted_params=cml.vertical.adj@coef,
                                     difference=my.data.combos.i$difference[k],difference_base=my.data.combos.i$difference_base[k])
      u.i <- u.i+u.ik
    }
    u_sq_sum <- u_sq_sum + u.i %*% t(u.i)
  }


  my_H <- - triangledown_u_sum/N
  my_J <- u_sq_sum/N
  my_G <- my_H %*% solve(my_J) %*% my_H
  my_variance_matrix <- solve(my_G)

  return(c(CLWPA_est=cml.vertical.adj@coef["theta"],
           CLWPA_se=((1/N)*my_variance_matrix[1,1])^.5,
           CLWPA_pval=2 * pnorm(abs((cml.vertical.adj@coef["theta"])/(((1/N)*my_variance_matrix[1,1])^.5)),
                                lower.tail = FALSE)))
}
