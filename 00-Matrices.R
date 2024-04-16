## Load packages
require(dplyr)

## Create the matrix A that converts between observations and DID estimators
### Several helper functions are used to create the components of A

gen_Aj <- function(j) {
  return(cbind(matrix(data=rep(-1, j-1), nrow=j-1, ncol=1),
               diag(nrow=j-1)))
}

gen_Adot <- function(J) {
  do.call(rbind,
          lapply(X=J:2,
                 FUN=function(x) cbind(matrix(data=0, nrow=x-1, ncol=J-x),
                                       gen_Aj(x))))
}

gen_Arow <- function(N,J,n1,n2,J2,Adot) {
  cbind(matrix(data=0, nrow=J2, ncol=(n1-1)*J),
        Adot,
        matrix(data=0, nrow=J2, ncol=(n2-n1-1)*J),
        -1*Adot,
        matrix(data=0, nrow=J2, ncol=(N-n2)*J))
}

gen_A <- function(N,J) {
  J2 <- J*(J-1)/2
  Adot <- gen_Adot(J)
  
  do.call(rbind,
          lapply(X=1:(N-1),
                 FUN=function(n1) do.call(rbind,
                                          lapply(X=(n1+1):N,
                                                 FUN=function(n2) gen_Arow(N,J,n1,n2,J2,Adot)))))
}

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

## Old versions of the gen_A function:
# gen_A <- function(N,J) {
#   J2 <- J*(J-1)/2
#   Adot <- gen_Adot(J)
#   
#   A <- NULL
#   for (n1 in 1:(N-1)) {
#     for (n2 in (n1+1):N) {
#       A <- rbind(A, gen_Arow(N,J,n1,n2,J2,Adot))
#     }
#   }
#   return(A)
# }
# 
# gen_A2 <- function(N,J) {
#   J2 <- J*(J-1)/2
#   Adot <- gen_Adot(J)
#   
#   A <- NULL
#   for (n1 in 1:(N-1)) {
#     A <- rbind(A, do.call(rbind,
#                           lapply(X=(n1+1):N,
#                                  FUN=function(n2) gen_Arow(N,J,n1,n2,J2,Adot))))
#   }
#   return(A)
# }

