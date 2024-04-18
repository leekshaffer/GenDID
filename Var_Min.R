#######################################
###### File: Var_Min.R ################
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/18 #############
###### Updated 2024/04/18 #############
#######################################

require(dplyr)

min_var <- function(DFT_obj,A_mat,v,solve_obj,Sigma) {
  if (nrow(Sigma) != ncol(Sigma) | nrow(Sigma) != ncol(A_mat)) {
    stop(simpleError("Sigma must be a square matrix with the same number of columns as A."))
  }
  base.w <- solve_obj$DID.weights %>% dplyr::select(starts_with("w.base"))
  Add.Obs.w <- solve_obj$DID.weights %>% dplyr::select(starts_with("Add.Obs.weights"))
  if (ncol(Add.Obs.w) == 0) {
    print("There is a unique solution to the constraint.")
    return(list(DID.weights=base.w,
                Obs.weights=solve_obj$Obs.weights %>% dplyr::select(starts_with("ATw.base"))))
  } else {
    dim_C <- ncol(Add.Obs.w)
    if (is.vector(v)) {
      v <- matrix(data=v, ncol=1)
      W_mat <- as.matrix(cbind(base.w,Add.Obs.w))
      Inn_mat <- t(W_mat) %*% A_mat %*% Sigma %*% t(A_mat) %*% W_mat
      opt_fn <- function(x) {matrix(data=c(1,x), nrow=1) %*% Inn_mat %*% matrix(data=c(1,x), ncol=1)}
      if (dim_C==1) { ## Solved quadratic case for k=1:
        if (Inn_mat[2,2] > 0) {
          Min <- -1*Inn_mat[1,2]/Inn_mat[2,2]
          Var <- as.numeric(opt_fn(Min))
        } else {
          stop(simpleError("No Minimum Found."))
        }
      } else {
        Opt <- optim(par=rep(0,dim_C),
                     fn=opt_fn,
                     method="Nelder-Mead")
        Min <- Opt$par
        Var <- Opt$value
      }
    } else {
      ## Need to iterate through each v
    }
  }
}