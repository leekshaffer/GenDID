#######################################
###### File: A_Const.R ################
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/16 #############
###### Updated 2024/04/16 #############
#######################################

## Load packages
require(dplyr)

## Create a data frame where each row corresponds to an entry in vector D
### The columns give i, i', j, and j' corresponding to that D
### Can be run either with integers N and J that represent the number of unique
#### clusters and periods, resp.
### Or with CL and PD, vectors giving the names/labels of the clusters and periods
### Or with a combination of CL or N and PD or J
gen_D <- function(CL=NULL,PD=NULL,N=NULL,J=NULL) {
  if (is.null(CL)) {
    if (is.null(N)) {
      stop(simpleError(message="Must specify either CL or N."))
    } else {
      if (is.numeric(N)) {
        if (N < 2) {
          stop(simpleError(message="N must be at least 2."))
        } else if (abs(N-round(N)) > .Machine$double.eps^2) {
          N <- round(N)
          warning(simpleWarning(message=paste0("N must be an integer. Using nearest integer: ",N,".")))
          CLu <- 1:N
        } else {
          CLu <- 1:N
        }
      } else {
        stop(simpleError(message="N must be a positive integer."))
      }
    }
  } else if (is.null(N)) {
    CLu <- unique(CL)
    N <- length(CLu)
    if (N < 2) {
      stop(simpleError(message="Must have at least 2 unique clusters."))
    }
  } else {
    CLu <- unique(CL)
    if (!is.null(N)) {
      if(is.numeric(N)) {
        if (N != length(CLu)) {
          warning(simpleWarning(message="N does not match number of unique clusters in CL. Using CL."))
        }
      }
    }
    N <- length(CLu)
    if (N < 2) {
      stop(simpleError(message="Must have at least 2 unique clusters."))
    }
  }
  
  if (is.null(PD)) {
    if (is.null(J)) {
      stop(simpleError(message="Must specify either PD or J."))
    } else {
      if (is.numeric(J)) {
        if (J < 2) {
          stop(simpleError(message="J must be at least 2."))
        } else if (abs(J-round(J)) > .Machine$double.eps^2) {
          J <- round(J)
          warning(simpleWarning(message=paste0("J must be an integer. Using nearest integer: ",J,".")))
          PDu <- 1:J
        } else {
          PDu <- 1:J
        }
      } else {
        stop(simpleError(message="J must be a positive integer."))
      }
    }
  } else if (is.null(J)) {
    PDu <- unique(PD)
    J <- length(PDu)
    if (J < 2) {
      stop(simpleError(message="Must have at least 2 unique time periods."))
    }
  } else {
    PDu <- unique(PD)
    if (!is.null(J)) {
      if(is.numeric(J)) {
        if (J != length(PDu)) {
          warning(simpleWarning(message="J does not match number of unique periods in PD. Using PD."))
        }
      }
    }
    J <- length(PDu)
    if (J < 2) {
      stop(simpleError(message="Must have at least 2 unique time periods."))
    }
  }
  
  Cl.df <- data.frame(i=CLu[unlist(sapply(1:(N-1), FUN=function(x) rep(x,N-x)))],
                      i.prime=CLu[unlist(sapply(1:(N-1), FUN=function(x) (x+1):N))])
  Pd.df <- data.frame(j=PDu[unlist(sapply(1:(J-1), FUN=function(x) rep(x,J-x)))],
                      j.prime=PDu[unlist(sapply(1:(J-1), FUN=function(x) (x+1):J))])
  return(list(Dmat=dplyr::cross_join(Cl.df, Pd.df),
              Dlen=N*(N-1)*J*(J-1)/4))
}
