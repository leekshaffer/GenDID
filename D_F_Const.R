#######################################
###### File: D_F_Const.R ##############
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/16 #############
###### Updated 2024/08/12 #############
#######################################

## Load packages
require(dplyr)
require(tibble)
require(tidyr)

## Create a data frame where each row corresponds to an entry in vector D
### The columns give i, i', j, and j' corresponding to that D
### Inputs are N, the number of clusters, and J, the number of periods.
gen_D <- function(N,J) {
  if (!is.numeric(N) | !is.numeric(J)) {
    stop(simpleError(message="N and J must both be numeric."))
  } else {
    if (N < 2 | J < 2) {
      stop(simpleError(message="N and J must both be at least 2."))
    } else {
      if (abs(N-round(N)) > .Machine$double.eps^2) {
        N <- round(N)
        warning(simpleWarning(message=paste0("N must be an integer. Using nearest integer: ",N,".")))
      }
      if (abs(J-round(J)) > .Machine$double.eps^2) {
        J <- round(J)
        warning(simpleWarning(message=paste0("J must be an integer. Using nearest integer: ",J,".")))
      }
      Cl.df <- data.frame(i=unlist(sapply(1:(N-1), FUN=function(x) rep(x,N-x))),
                          i.prime=unlist(sapply(1:(N-1), FUN=function(x) (x+1):N)))
      Pd.df <- data.frame(j=unlist(sapply(1:(J-1), FUN=function(x) rep(x,J-x))),
                          j.prime=unlist(sapply(1:(J-1), FUN=function(x) (x+1):J)))
      return(dplyr::cross_join(Cl.df, Pd.df))
    }
  }
}

## Take a vector of unique clusters and a corresponding vector of their first periods on intervention
### Note that Clusters and StartPeriods can be labelled however desired
### OrderedPds is a vector of periods where outcomes will be recorded, in order
gen_js <- function(Clusters,StartPeriods,OrderedPds) {
  if (length(unique(OrderedPds)) != length(OrderedPds)) {
    warning(simpleWarning("There are repeated names in OrderedPds. The first occurrence only will be used."))
    OrderedPds <- unique(OrderedPds)
  }
  J <- length(OrderedPds)
  
  if (length(Clusters) != length(unique(Clusters))) {
    stop(simpleError("Each cluster must appear only once."))
  }
  N <- length(unique(Clusters))
  if (N < 2) {
    stop(simpleError("There must be at least 2 clusters."))
  } else if (N != length(StartPeriods)) {
    stop(simpleError("The length of Clusters and StartPeriods must match."))
  }
  
  PdLabels <- tibble(Labels=OrderedPds, Pd.Num=1:J)
  Start_js <- tibble(Clusters=Clusters,
                     Labels=StartPeriods)
  if (sum(StartPeriods %in% OrderedPds) < length(StartPeriods)) {
    if (is.numeric(OrderedPds) & is.numeric(StartPeriods) &
        sum(OrderedPds[order(OrderedPds)]==OrderedPds)==J) {
      warning(simpleWarning("Since OrderedPds is numeric, StartPeriods outside this range are treated numerically. NAs/Infs are considered never-treated."))
      min_Ord <- min(OrderedPds)
      Start_js <- Start_js %>% 
        dplyr::mutate(Pd.Num=if_else(is.na(Labels) | is.infinite(Labels), 
                                     Inf,
                                     Labels-min_Ord+1))
    } else {
      warning(simpleWarning("All StartPeriods not in OrderedPds are considered never-treated units within this setting."))
      Start_js <- Start_js %>%
        left_join(PdLabels, by=join_by(Labels)) %>%
        dplyr::mutate(Pd.Num=if_else(is.na(Pd.Num),Inf,Pd.Num))
    }
  } else {
    Start_js <- Start_js %>%
      dplyr::left_join(PdLabels, by=join_by(Labels))
  }
  
  # Start_js <- tibble(Clusters=Clusters,
  #                    Labels=StartPeriods) %>%
  #   left_join(PdLabels, by=join_by(Labels)) %>%
  #   dplyr::mutate(Pd.Num=if_else(is.na(Pd.Num),Inf,Pd.Num)) %>%
  Start_js <- Start_js %>%
    dplyr::arrange(Pd.Num) %>%
    dplyr::mutate(Cl.Num=1:N) %>%
    dplyr::rename(Start_j = Pd.Num) %>%
    dplyr::select(Clusters,Cl.Num,Start_j)
  
  return(list(Start_js=Start_js,
              PdLabels=PdLabels %>% dplyr::rename(Periods=Labels) %>% 
                dplyr::select(Pd.Num,Periods),
              N=N,
              J=J))
}

gen_Theta <- function(js_Obj,Assumption) {
  Theta <- js_Obj$Start_js %>% 
    dplyr::cross_join(js_Obj$PdLabels) %>%
    dplyr::filter(Pd.Num >= Start_j) %>%
    mutate(Diff=Pd.Num - Start_j + 1)
  
  if (Assumption==5) {
    JoinBy <- NULL
  } else if (Assumption==4) {
    JoinBy <- c("Pd.Num")
  } else if (Assumption==3) {
    JoinBy <- c("Diff")
  } else if (Assumption==2) {
    JoinBy=c("Pd.Num","Diff")
  } else if (Assumption==1) {
    JoinBy=c("Cl.Num","Pd.Num","Diff")
  } else {
    stop(simpleError("Assumption must be a value 1 through 5 corresponding to the assumption setting desired."))
  }
  
  if(is.null(JoinBy)) {
    All <- tibble(Theta=1)
  } else {
    All <- Theta %>% dplyr::select(all_of(JoinBy)) %>% 
      dplyr::arrange(across(JoinBy)) %>%
      dplyr::distinct() %>%
      mutate(Theta=row_number())
  }
  Full <- merge(Theta, All, by=JoinBy)
  
  Schematic <- matrix(data=0,
                      nrow=js_Obj$N,
                      ncol=js_Obj$J)
  Schematic[as.matrix(Full %>% dplyr::select(Cl.Num,Pd.Num))] <- Full$Theta
  dimnames(Schematic) <- list(js_Obj$Start_js$Clusters,js_Obj$PdLabels$Periods)
  
  return(list(All=All,
         Full=Full,
         Schematic=Schematic))
}

gen_F <- function(D,Theta) {
  ThetaF <- Theta$Full %>% dplyr::select(Cl.Num,Pd.Num,Theta)
  D_aug <- D %>% dplyr::select(i,i.prime,j,j.prime,Type,TypeLabel) %>%
    dplyr::left_join(ThetaF %>% 
                       dplyr::rename(Pos.i=Theta),
                     by=join_by(i==Cl.Num,
                                j.prime==Pd.Num)) %>%
    dplyr::left_join(ThetaF %>%
                       dplyr::rename(Neg.i=Theta),
                     by=join_by(i==Cl.Num,
                                j==Pd.Num)) %>%
    dplyr::left_join(ThetaF %>% 
                       dplyr::rename(Neg.i.prime=Theta),
                     by=join_by(i.prime==Cl.Num,
                                j.prime==Pd.Num)) %>%
    dplyr::left_join(ThetaF %>%
                       dplyr::rename(Pos.i.prime=Theta),
                     by=join_by(i.prime==Cl.Num,
                                j==Pd.Num)) %>%
    mutate(Row=row_number())
  
  F_mat <- matrix(0,nrow=dim(D_aug)[1], ncol=length(Theta$All$Theta))
  F_mat[as.matrix(D_aug %>% dplyr::select(Row,Pos.i) %>% dplyr::filter(!is.na(Pos.i)))] <- 
    F_mat[as.matrix(D_aug %>% dplyr::select(Row,Pos.i) %>% dplyr::filter(!is.na(Pos.i)))]+1
  F_mat[as.matrix(D_aug %>% dplyr::select(Row,Neg.i) %>% dplyr::filter(!is.na(Neg.i)))] <- 
    F_mat[as.matrix(D_aug %>% dplyr::select(Row,Neg.i) %>% dplyr::filter(!is.na(Neg.i)))]-1
  F_mat[as.matrix(D_aug %>% dplyr::select(Row,Neg.i.prime) %>% dplyr::filter(!is.na(Neg.i.prime)))] <- 
    F_mat[as.matrix(D_aug %>% dplyr::select(Row,Neg.i.prime) %>% dplyr::filter(!is.na(Neg.i.prime)))]-1
  F_mat[as.matrix(D_aug %>% dplyr::select(Row,Pos.i.prime) %>% dplyr::filter(!is.na(Pos.i.prime)))] <- 
    F_mat[as.matrix(D_aug %>% dplyr::select(Row,Pos.i.prime) %>% dplyr::filter(!is.na(Pos.i.prime)))]+1
  return(list(D_aug=D_aug, F_mat=F_mat))
}

gen_DFT <- function(Clusters, StartPeriods, OrderedPds,
                    Assumption=5) {
  js_Obj <- gen_js(Clusters, StartPeriods, OrderedPds)
  
  Theta <- gen_Theta(js_Obj, Assumption)
  
  D <- gen_D(js_Obj$N,js_Obj$J) %>% 
    left_join(js_Obj$Start_js %>%
                dplyr::rename(i=Cl.Num,Start.i=Start_j,
                              Cl.i=Clusters), by="i") %>%
    left_join(js_Obj$Start_js %>% 
                dplyr::rename(i.prime=Cl.Num,Start.i.prime=Start_j,
                              Cl.i.prime=Clusters), by="i.prime") %>%
    mutate(i.j.leadlag=j-Start.i+1,
           ip.j.leadlag=j-Start.i.prime+1,
           i.jp.leadlag=j.prime-Start.i+1,
           ip.jp.leadlag=j.prime-Start.i.prime+1,
           i.type=if_else(i.j.leadlag > 0,"Always-Treated",
                          if_else(i.jp.leadlag <= 0, "Always-Control","Switch")),
           ip.type=if_else(ip.j.leadlag > 0,"Always-Treated",
                          if_else(ip.jp.leadlag <= 0, "Always-Control","Switch")),
           Type=if_else(i.type=="Always-Control",1,
                        if_else(i.type=="Switch",
                                if_else(ip.type=="Switch",4,2),
                                if_else(ip.type=="Always-Control",3,
                                        if_else(ip.type=="Switch",5,6)))),
           TypeLabel=factor(Type,
           levels=1:6, labels=c("Both Always-Control",
                                "Switch vs. Always-Control",
                                "Always-Treated vs. Always-Control",
                                "Both Switch",
                                "Always-Treated vs. Switch",
                                "Both Always-Treated")))
  
  Fres <- gen_F(D, Theta)
  return(list(D_aug=Fres$D_aug, F_mat=Fres$F_mat, Theta=Theta, 
              N=js_Obj$N, J=js_Obj$J))
}