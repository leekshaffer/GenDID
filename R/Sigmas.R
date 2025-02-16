#######################################
###### File: Sigmas.R #################
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/19 #############
###### Updated 2024/04/19 #############
#######################################

## Possibilities to speed up diag matrix creation:
#### Matrix pkg, bdiag function (sparse matrix?)




#' Convert Cluster Covariance to Full Block-Diagonal Matrix
#'
#' This function takes a covariance matrix for a single cluster (\eqn{\Sigma_i}) and
#' constructs a full block-diagonal covariance matrix for \eqn{N} exchangeable clusters.
#'
#' @param Sigma.i A square matrix representing the covariance structure of a single cluster.
#' @param N An integer specifying the number of clusters.
#'
#' @return A block-diagonal matrix of dimensions \eqn{(N \times J) \times (N \times J)},
#' where \eqn{J} is the number of rows (or columns) in \eqn{\Sigma_i}.
#'
#' @details
#' This function assumes that all clusters are exchangeable, meaning they share the
#' same covariance structure. It creates a block-diagonal matrix where each block
#' is the input \eqn{\Sigma_i} repeated \eqn{N} times.
#'
#' @examples
#' # Example usage:
#' Sigma.i <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' N <- 3
#' full_matrix <- cluster_to_full(Sigma.i, N)
#' print(full_matrix)
#'
#' @export
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




#' Convert Correlation Matrix to Covariance Matrix
#'
#' This function takes a correlation matrix and a vector of standard deviations
#' to generate a covariance matrix.
#'
#' @param Corr A square matrix representing the correlation structure.
#' @param SDvec A numeric vector of standard deviations. Its length must either
#' match the dimensions of \code{Corr}, be 1 (a constant value for all variables),
#' or be repeated/truncated as needed.
#'
#' @return A covariance matrix of the same dimensions as the input correlation matrix.
#' If \code{SDvec} is \code{NULL}, the input correlation matrix is returned unchanged.
#'
#' @details
#' The function computes the covariance matrix using the formula:
#' \deqn{\Sigma = \text{diag(SDvec)} \cdot \text{Corr} \cdot \text{diag(SDvec)}}
#' where \code{diag(SDvec)} creates a diagonal matrix with elements of \code{SDvec}.
#'
#' If the length of \code{SDvec} is 1, all variables are assumed to have the same
#' standard deviation. If the length of \code{SDvec} does not match the dimensions
#' of \code{Corr}, it will be repeated or truncated to fit, with a warning issued.
#'
#' @examples
#' # Example 1: Correlation matrix with identical SDs
#' Corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' SDvec <- 2
#' cov_matrix <- corr_to_cov(Corr, SDvec)
#' print(cov_matrix)
#'
#' # Example 2: Correlation matrix with different SDs
#' SDvec <- c(2, 3)
#' cov_matrix <- corr_to_cov(Corr, SDvec)
#' print(cov_matrix)
#'
#' # Example 3: Using the correlation matrix directly
#' cov_matrix <- corr_to_cov(Corr, NULL)
#' print(cov_matrix)
#'
#' @export
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



#' Create Covariance Matrix for Independent Time Periods
#'
#' This function generates a covariance matrix for independent time periods, given the number of clusters (\code{N}),
#' the number of time periods per cluster (\code{J}), and optionally a vector of standard deviations (\code{SDvec}).
#'
#' @param N An integer representing the number of clusters.
#' @param J An integer representing the number of time periods within each cluster.
#' @param SDvec An optional numeric vector of standard deviations. If \code{NULL}, all variances are assumed to be 1.
#'
#' @return A covariance matrix of dimensions \eqn{(N \times J) \times (N \times J)}:
#' \describe{
#'   \item{If \code{SDvec} is \code{NULL}}{The function returns an identity matrix, assuming unit variance for all periods.}
#'   \item{If \code{length(SDvec) == J}}{The function returns a block-diagonal covariance matrix, where each block corresponds to \code{diag(SDvec^2)} repeated across clusters.}
#'   \item{Otherwise}{The function assumes the time periods are independent and calls \code{cluster_to_full} to construct the covariance matrix based on the given \code{SDvec}.}
#' }
#'
#' @details
#' This function handles scenarios where:
#' \itemize{
#'   \item Variances are uniform across all clusters and time periods.
#'   \item Variances are specific to each time period but consistent across clusters.
#'   \item Variances are both cluster- and time-specific.
#' }
#' The function uses \code{diag} to create diagonal matrices and leverages the helper functions \code{cluster_to_full} and \code{corr_to_cov} as needed.
#'
#' @examples
#' # Example 1: Identity covariance matrix for 2 clusters and 3 time periods
#' Sigma <- create_Sigma_Ind(N = 2, J = 3)
#' print(Sigma)
#'
#' # Example 2: Covariance matrix with specified SDvec for 2 clusters and 3 time periods
#' SDvec <- c(1, 2, 3)
#' Sigma <- create_Sigma_Ind(N = 2, J = 3, SDvec = SDvec)
#' print(Sigma)
#'
#' @export
create_Sigma_Ind <- function(N,J,SDvec=NULL) {
  if (is.null(SDvec)) {
    return(diag(1, nrow=N*J))
  } else if (length(SDvec)==J) {
    return(diag(x=SDvec^2, nrow=N*J))
  } else {
    return(cluster_to_full(corr_to_cov(diag(1,nrow=J),SDvec),N))
  }
}

#' Create Covariance Matrix for Exchangeable Time Periods
#'
#' This function generates a covariance matrix for exchangeable time periods (compound symmetry structure)
#' based on the correlation between time periods (\code{rho}), the number of clusters (\code{N}),
#' the number of time periods per cluster (\code{J}), and an optional vector of standard deviations (\code{SDvec}).
#'
#' @param rho A numeric value representing the correlation between time periods.
#' @param N An integer representing the number of clusters.
#' @param J An integer representing the number of time periods within each cluster.
#' @param SDvec An optional numeric vector of standard deviations. If \code{NULL}, unit variances are assumed.
#'
#' @return A covariance matrix of dimensions \eqn{(N \times J) \times (N \times J)}:
#' \describe{
#'   \item{Correlation matrix (\eqn{\rho})}{The off-diagonal elements of the covariance structure are set to \eqn{\rho}.}
#'   \item{Diagonal elements}{Variances are derived from \code{SDvec} if provided; otherwise, unit variances are used.}
#'   \item{Block structure}{The matrix is constructed as a block-diagonal matrix, repeating the same structure for each cluster.}
#' }
#'
#' @details
#' This function assumes a compound symmetry structure for time periods within each cluster, where:
#' \itemize{
#'   \item All off-diagonal correlations between time periods are equal to \code{rho}.
#'   \item Variances can be uniform or specified individually using \code{SDvec}.
#' }
#' The function leverages \code{cluster_to_full} and \code{corr_to_cov} to construct the covariance matrix.
#'
#' @examples
#' # Example 1: Exchangeable time periods with rho = 0.5, 2 clusters, and 3 periods
#' Sigma <- create_Sigma_CS(rho = 0.5, N = 2, J = 3)
#' print(Sigma)
#'
#' # Example 2: Exchangeable time periods with rho = 0.3 and specified SDvec
#' SDvec <- c(1, 2, 3)
#' Sigma <- create_Sigma_CS(rho = 0.3, N = 2, J = 3, SDvec = SDvec)
#' print(Sigma)
#'
#' @export
create_Sigma_CS <- function(rho, N, J, SDvec=NULL) {
  return(cluster_to_full(
    corr_to_cov(matrix(rho, nrow=J, ncol=J) + diag(1-rho, nrow=J),
                SDvec), N))
}

#' Create Covariance Matrix for Auto-Regressive Time Periods
#'
#' This function generates a covariance matrix for auto-regressive time periods of order 1 (AR(1)),
#' based on the correlation decay parameter (\code{rho}), the number of clusters (\code{N}),
#' the number of time periods per cluster (\code{J}), and an optional vector of standard deviations (\code{SDvec}).
#'
#' @param rho A numeric value (0 ≤ \code{rho} ≤ 1) representing the auto-regressive correlation between adjacent time periods.
#' @param N An integer specifying the number of clusters.
#' @param J An integer specifying the number of time periods within each cluster.
#' @param SDvec An optional numeric vector of standard deviations. If \code{NULL}, unit variances are assumed.
#'
#' @return A covariance matrix of dimensions \eqn{(N \times J) \times (N \times J)}:
#' \describe{
#'   \item{Correlation structure}{The correlation between time periods decays geometrically based on their separation, \eqn{\rho^{|t-s|}}.}
#'   \item{Block structure}{The matrix is constructed as a block-diagonal matrix, repeating the same structure for each cluster.}
#'   \item{Variance scaling}{If \code{SDvec} is provided, variances are adjusted accordingly.}
#' }
#'
#' @details
#' This function assumes an AR(1) structure for time periods within each cluster, where:
#' \itemize{
#'   \item The correlation between two time periods \eqn{t} and \eqn{s} is \eqn{\rho^{|t-s|}}.
#'   \item Variances can be uniform or specified using \code{SDvec}.
#' }
#' The function constructs the AR(1) correlation matrix for one cluster and then uses \code{cluster_to_full}
#' and \code{corr_to_cov} to scale it with variances and replicate it across clusters.
#'
#' @examples
#' # Example 1: AR(1) structure with rho = 0.5, 2 clusters, and 3 periods
#' Sigma <- create_Sigma_AR1(rho = 0.5, N = 2, J = 3)
#' print(Sigma)
#'
#' # Example 2: AR(1) structure with rho = 0.7 and specified SDvec
#' SDvec <- c(1, 2, 3)
#' Sigma <- create_Sigma_AR1(rho = 0.7, N = 2, J = 3, SDvec = SDvec)
#' print(Sigma)
#'
#' @export
create_Sigma_AR1 <- function(rho,N,J,SDvec=NULL) {
  Corr_Cl <- do.call(rbind, lapply(1:J, FUN=function(i) rho^(abs((1-i):(J-i)))))
  return(cluster_to_full(corr_to_cov(Corr_Cl,SDvec), N))
}
