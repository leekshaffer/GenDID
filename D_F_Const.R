#######################################
###### File: A_Const.R ################
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/16 #############
###### Updated 2024/04/18 #############
#######################################

## Load packages
require(dplyr)
require(tibble)
require(tidyr)

## Create a data frame where each row corresponds to an entry in vector D
### The columns give i, i', j, and j' corresponding to that D
### Can be run either with integers N and J that represent the number of unique
#### clusters and periods, resp.
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

## Take a vector of unique clusters and a corresponding vector of their first periods on intervention,
### Plus output from gen_D and return types for E[D]
### Note that Clusters and StartPeriods can be labelled however desired
### If list of periods is anything other than 1:J, 
#### a vector giving the ordered period names should be included as PeriodOrder
### Otherwise J should indicate the largest period
gen_Start_js <- function(Clusters,StartPeriods,J=NULL,PeriodOrder=NULL) {
  N <- length(Clusters)
  if (N < 2) {
    stop(simpleError("There must be at least 2 clusters."))
  } else if (N != length(StartPeriods)) {
    stop(simpleError("The length of Clusters and StartPeriods must match."))
  }
  
  if (!is.null(PeriodOrder)) {
    if (length(unique(PeriodOrder) != length(PeriodOrder))) {
      warning(simpleWarning("There are repeated names in PeriodOrder. The first occurrence only will be used."))
      PeriodOrder <- unique(PeriodOrder)
    }
    if (!is.null(J)) {
      if (J != length(PeriodOrder)) {
        warning(simpleWarning("J is not the total number of periods, ignoring J"))
        J <- length(PeriodOrder)
      }
    } else {
      J <- length(PeriodOrder)
    }
    OrderedPds <- tibble(Labels=PeriodOrder, Periods=1:J)
    Start_js <- tibble(Clusters=Clusters,Labels=StartPeriods) %>% 
      left_join(OrderedPds,by="Labels") %>% 
      dplyr::rename(Start_j=Periods) %>%
      dplyr::select(c(Clusters,Start_j)) %>%
      dplyr::arrange(Start_j) %>%
      mutate(Cl.Num=1:N)
  } else if (is.null(J)) {
    stop(simpleError("Must specify either the vector of ordered period names PeriodOrder or the total number of periods J"))
  } else if (!is.numeric(J)) {
    stop(simpleError("J must be numeric"))
  } else if (J < 2) {
    stop(simpleError("J must be at least 2."))
  } else {
    if (abs(J-round(J)) > .Machine$double.eps^2) {
      J <- round(J)
      warning(simpleWarning(paste0("J must be an integer. Using nearest integer: ",J)))
    }
    if (!is.numeric(StartPeriods)) {
      stop(simpleError("Must either specify numeric StartPeriods or ordered PeriodOrder."))
    } else if (min(StartPeriods) < 1 | max(StartPeriods) > J) {
      stop(simpleError("All StartPeriods must be between 1 and J, inclusive."))
    } else {
      Start_js <- tibble(Clusters=Clusters,Start_j=StartPeriods) %>%
        dplyr::arrange(Start_j) %>%
        mutate(Cl.Num=1:N)
      OrderedPds <- tibble(Labels=1:J, Periods=1:J)
    }
  }
  return(list(Start_js=Start_js,OrderedPds=OrderedPds))
}

gen_Theta <- function(Start_js,OrderedPds,Assumption) {
  Theta <- Start_js %>% dplyr::cross_join(OrderedPds) %>%
    dplyr::filter(Periods >= Start_j) %>%
    mutate(Diff=Periods - Start_j + 1)
  
  if (Assumption==5) {
    JoinBy <- NULL
  } else if (Assumption==4) {
    JoinBy <- c("Periods")
  } else if (Assumption==3) {
    JoinBy <- c("Diff")
  } else if (Assumption==2) {
    JoinBy=c("Periods","Diff")
  } else if (Assumption==1) {
    JoinBy=c("Cl.Num","Periods","Diff")
  } else {
    stop(simpleError("Assumption must be a value 1 through 5 corresponding to the assumption setting desired."))
  }
  
  if(is.null(JoinBy)) {
    All <- tibble(Theta=1)
    Full <- Theta %>%
      dplyr::cross_join(All)
  } else {
    All <- Theta %>% dplyr::select(all_of(JoinBy)) %>% 
      dplyr::arrange(across(JoinBy)) %>%
      dplyr::distinct() %>%
      mutate(Theta=row_number())
    Full <- Theta %>% 
      dplyr::left_join(All, by=JoinBy)
  }
  
  Schematic <- matrix(data=0,
                      nrow=length(Start_js$Cl.Num),
                      ncol=length(OrderedPds$Periods))
  Schematic[as.matrix(Full %>% dplyr::select(Cl.Num,Periods))] <- Full$Theta
  dimnames(Schematic) <- list(Start_js$Clusters,OrderedPds$Labels)
  
  return(list(All=All,
         Full=Full,
         Schematic=Schematic))
}

gen_F <- function(D,Theta) {
  ThetaF <- Theta$Full %>% dplyr::select(Cl.Num,Periods,Theta)
  D_aug <- D %>% dplyr::select(i,i.prime,j,j.prime,Type,TypeLabel) %>%
    dplyr::left_join(ThetaF %>% 
                       dplyr::rename(i=Cl.Num,
                                     j.prime=Periods,
                                     Pos.i=Theta)) %>%
    dplyr::left_join(ThetaF %>%
                       dplyr::rename(i=Cl.Num,
                                     j=Periods,
                                     Neg.i=Theta)) %>%
    dplyr::left_join(ThetaF %>% 
                       dplyr::rename(i.prime=Cl.Num,
                                     j.prime=Periods,
                                     Neg.i.prime=Theta)) %>%
    dplyr::left_join(ThetaF %>%
                       dplyr::rename(i.prime=Cl.Num,
                                     j=Periods,
                                     Pos.i.prime=Theta)) %>%
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

gen_DFT <- function(Clusters,StartPeriods,J=NULL,PeriodOrder=NULL,Assumption=5) {
  SJ <- gen_Start_js(Clusters, StartPeriods, J, PeriodOrder)
  
  Theta <- gen_Theta(SJ$Start_js, SJ$OrderedPds, Assumption)
  
  D <- gen_D(length(SJ$Start_js$Clusters),J) %>% 
    left_join(SJ$Start_js %>%
                dplyr::rename(i=Cl.Num,Start.i=Start_j,
                              Cl.i=Clusters), by="i") %>%
    left_join(SJ$Start_js %>% 
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
              N=length(SJ$Start_js$Clusters),J=J))
}



### Old gen_D function:

# gen_D <- function(CL=NULL,PD=NULL,N=NULL,J=NULL) {
#   if (is.null(CL)) {
#     if (is.null(N)) {
#       stop(simpleError(message="Must specify either CL or N."))
#     } else {
#       if (is.numeric(N)) {
#         if (N < 2) {
#           stop(simpleError(message="N must be at least 2."))
#         } else if (abs(N-round(N)) > .Machine$double.eps^2) {
#           N <- round(N)
#           warning(simpleWarning(message=paste0("N must be an integer. Using nearest integer: ",N,".")))
#           CLu <- 1:N
#         } else {
#           CLu <- 1:N
#         }
#       } else {
#         stop(simpleError(message="N must be a positive integer."))
#       }
#     }
#   } else if (is.null(N)) {
#     CLu <- unique(CL)
#     N <- length(CLu)
#     if (N < 2) {
#       stop(simpleError(message="Must have at least 2 unique clusters."))
#     }
#   } else {
#     CLu <- unique(CL)
#     if (!is.null(N)) {
#       if(is.numeric(N)) {
#         if (N != length(CLu)) {
#           warning(simpleWarning(message="N does not match number of unique clusters in CL. Using CL."))
#         }
#       }
#     }
#     N <- length(CLu)
#     if (N < 2) {
#       stop(simpleError(message="Must have at least 2 unique clusters."))
#     }
#   }
#   
#   if (is.null(PD)) {
#     if (is.null(J)) {
#       stop(simpleError(message="Must specify either PD or J."))
#     } else {
#       if (is.numeric(J)) {
#         if (J < 2) {
#           stop(simpleError(message="J must be at least 2."))
#         } else if (abs(J-round(J)) > .Machine$double.eps^2) {
#           J <- round(J)
#           warning(simpleWarning(message=paste0("J must be an integer. Using nearest integer: ",J,".")))
#           PDu <- 1:J
#         } else {
#           PDu <- 1:J
#         }
#       } else {
#         stop(simpleError(message="J must be a positive integer."))
#       }
#     }
#   } else if (is.null(J)) {
#     PDu <- unique(PD)
#     J <- length(PDu)
#     if (J < 2) {
#       stop(simpleError(message="Must have at least 2 unique time periods."))
#     }
#   } else {
#     PDu <- unique(PD)
#     if (!is.null(J)) {
#       if(is.numeric(J)) {
#         if (J != length(PDu)) {
#           warning(simpleWarning(message="J does not match number of unique periods in PD. Using PD."))
#         }
#       }
#     }
#     J <- length(PDu)
#     if (J < 2) {
#       stop(simpleError(message="Must have at least 2 unique time periods."))
#     }
#   }
#   
#   Cl.df <- data.frame(i=CLu[unlist(sapply(1:(N-1), FUN=function(x) rep(x,N-x)))],
#                       i.prime=CLu[unlist(sapply(1:(N-1), FUN=function(x) (x+1):N))])
#   Pd.df <- data.frame(j=PDu[unlist(sapply(1:(J-1), FUN=function(x) rep(x,J-x)))],
#                       j.prime=PDu[unlist(sapply(1:(J-1), FUN=function(x) (x+1):J))])
#   return(list(Dmat=dplyr::cross_join(Cl.df, Pd.df),
#               Dlen=N*(N-1)*J*(J-1)/4))
# }


