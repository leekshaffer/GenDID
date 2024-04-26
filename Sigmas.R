#######################################
###### File: Sigmas.R #################
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/19 #############
###### Updated 2024/04/19 #############
#######################################

## Possibilities to speed up diag matrix creation:
#### Matrix pkg, bdiag function (sparse matrix?)


#### Methods for exchangeable clusters ####
### Function that takes Sigma.i for one cluster and creates full block-diagonal
cluster_to_full <- function(Sigma.i,N) {
  if (nrow(Sigma.i) != ncol(Sigma.i)) {
    stop(simpleError("Sigma.i must be a square matrix"))
  }
  J <- nrow(Sigma.i)
  Sigma <- matrix(0, nrow=N*J, ncol=N*J)
  for (i in 1:N) {
    Sigma[((i-1)*J+1):(i*J),((i-1)*J+1):(i*J)] <- Sigma.i
  }
  return(Sigma)
}
### Function that takes a correlation matrix and a SDvec and creates covariance matrix
corr_to_cov <- function(Corr,SDvec) {
  if (nrow(Corr) != ncol(Corr)) {
    stop(simpleError("The correlation matrix must be square, with dimension equal to the length of the SD vector."))
  } else if (is.null(SDvec)) {
    return(Corr)
  } else if (length(SDvec)==1) {
    return(SDvec^2*Corr)
  } else if (length(unique(SDvec))==1) {
    return(SDvec[1]^2*Corr)
  } else if (length(SDvec)==nrow(Corr)) {
    return(diag(x=SDvec) %*% Corr %*% diag(SDvec))
  } else if (length(SDvec) != nrow(Corr)) {
    warning(simpleWarning("The length of SDvec is not the same as the dimension of Corr. It will be repeated/truncated as necessary."))
    return(diag(x=rep(SDvec, length.out=nrow(Corr))) %*% Corr %*% diag(x=rep(SDvec, length.out=nrow(Corr))))
  } else {
    stop(simpleError("Incompatible arguments."))
  }
}

### Independent time periods ###
create_Sigma_Ind <- function(N,J,SDvec=NULL) {
  if (is.null(SDvec)) {
    return(diag(1, nrow=N*J))
  } else if (length(SDvec)==J) {
    return(diag(x=SDvec^2, nrow=N*J))
  } else {
    return(cluster_to_full(corr_to_cov(diag(1,nrow=J),SDvec),N))
  }
}

### Exchangeable time periods ###
create_Sigma_CS <- function(rho,N,J,SDvec=NULL) {
  return(cluster_to_full(
    corr_to_cov(matrix(rho, nrow=J, ncol=J) + diag(1-rho, nrow=J),
                SDvec),N))
}

### Auto-Regressive time periods ###
create_Sigma_AR1 <- function(rho,N,J,SDvec=NULL) {
  Corr_Cl <- do.call(rbind, lapply(1:J, FUN=function(i) rho^(abs((1-i):(J-i)))))
  return(cluster_to_full(corr_to_cov(Corr_Cl,SDvec), N))
}
