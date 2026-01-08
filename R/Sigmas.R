#######################################
###### File: Sigmas.R #################
#######################################

# cluster_to_full function

### Inputs:
#### Sigma.i: A square matrix representing the variance-covariance matrix of a single unit
#### N: An integer specifying the number of units
### Output:
#### A block-diagonal matrix with Sigma.i repeated N times on the diagonal, and 0s elsewhere

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

# corr_to_cov function

### Inputs:
#### Corr: A square matrix representing the correlation structure.
#### SDvec: A numeric vector of standard deviations
####  length must be the dimension of Corr or 1 (constant value for all variables)
####  else it will be repeated/truncated
### Output:
#### A variance-covariance matrix of the same dimensions as the input correlation matrix.

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

# create_Simga_Ind function

### Inputs:
#### N: An integer representing the number of units.
#### J: An integer representing the number of time periods.
#### SDvec: An optional numeric vector of standard deviations.
### Output:
#### A variance-covariance matrix of dimensions NJ by NJ, with independent time periods and the specified SDs.

create_Sigma_Ind <- function(N,J,SDvec=NULL) {
  if (is.null(SDvec)) {
    return(diag(1, nrow=N*J))
  } else if (length(SDvec)==J) {
    return(diag(x=SDvec^2, nrow=N*J))
  } else {
    return(cluster_to_full(corr_to_cov(diag(1,nrow=J),SDvec),N))
  }
}

# create_Simga_CS function

### Inputs:
#### rho: A numeric value representing the within-unit correlation between any two time periods.
#### N: An integer representing the number of units.
#### J: An integer representing the number of time periods.
#### SDvec: An optional numeric vector of standard deviations.
### Output:
#### A variance-covariance matrix of dimensions NJ by NJ,
####  with exchangeable (compound symmetric) within-unit correlations and the specified SDs.

create_Sigma_CS <- function(rho, N, J, SDvec=NULL) {
  return(cluster_to_full(
    corr_to_cov(matrix(rho, nrow=J, ncol=J) + diag(1-rho, nrow=J),
                SDvec), N))
}

# create_Simga_AR1 function

### Inputs:
#### rho: A numeric value representing the auto-regressive within-unit correlation
####  between two time periods that are one period apart.
#### N: An integer representing the number of units.
#### J: An integer representing the number of time periods.
#### SDvec: An optional numeric vector of standard deviations.
### Output:
#### A variance-covariance matrix of dimensions NJ by NJ,
####  with auto-regressive (lag 1) within-unit correlations and the specified SDs.

create_Sigma_AR1 <- function(rho,N,J,SDvec=NULL) {
  Corr_Cl <- do.call(rbind, lapply(1:J, FUN=function(i) rho^(abs((1-i):(J-i)))))
  return(cluster_to_full(corr_to_cov(Corr_Cl,SDvec), N))
}
