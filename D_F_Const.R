#######################################
###### File: A_Const.R ################
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/16 #############
###### Updated 2024/04/16 #############
#######################################

## Load packages
require(dplyr)
require(tibble)

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
gen_ED <- function(Clusters,StartPeriods,J=NULL,PeriodOrder=NULL) {
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
    }
  }
  
  D <- gen_D(N,J) %>% 
    left_join(Start_js %>%
                dplyr::rename(i=Cl.Num,Start.i=Start_j,
                              Cl.i=Clusters), by="i") %>%
    left_join(Start_js %>% 
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
           Type=factor(if_else(i.type=="Always-Control",1,
                        if_else(i.type=="Switch",
                                if_else(ip.type=="Switch",4,2),
                                if_else(ip.type=="Always-Control",3,
                                        if_else(ip.type=="Switch",5,6)))),
           levels=1:6, labels=c("Both Always-Control",
                                "Switch vs. Always-Control",
                                "Always-Treated vs. Always-Control",
                                "Both Switch",
                                "Always-Treated vs. Switch",
                                "Both Always-Treated")))
  
  Theta <- Start_js %>% dplyr::cross_join(tibble(Periods=1:J)) %>%
    dplyr::filter(Periods >= Start_j) %>%
    mutate(Diff=Periods - Start_j + 1)
  
  if (Assumption==5) {
    Theta_All <- Theta %>% mutate(Theta=1)
  } else if (Assumption==4) {
    Theta_All <- Theta %>% 
      dplyr::left_join(Theta %>% dplyr::select(Periods) %>% 
                         dplyr::distinct() %>%
                         mutate(Theta=row_number()),
                       by="Periods")
  } else if (Assumption==3) {
    Theta_All <- Theta %>% 
      dplyr::left_join(Theta %>% dplyr::select(Diff) %>% 
                         dplyr::distinct() %>%
                         mutate(Theta=row_number()),
                       by="Diff")
  } else if (Assumption==2) {
    Theta_All <- Theta %>%
      dplyr::left_join()
  } else if (Assumption==1) {
    
  } else {
    stop(simpleError("Assumption must be a value 1 through 5 corresponding to the assumption setting desired."))
  }
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


