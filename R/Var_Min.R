#######################################
###### File: Var_Min.R ################
###### Lee Kennedy-Shaffer ############
#######################################

# Helper functions

## opt_fn function

### Inputs:
#### x: the vector, which with a preceding 1, is on the outsides of the quadratic form
#### Inner: the square matrix which is inside the quadratic form, of dimension 1 more than x
### Output: matrix resulting from quadratic form multiplication

opt_fn <- function(x,Inner) {
  matrix(data=c(1,x), nrow=1) %*% Inner %*% matrix(data=c(1,x), ncol=1)
}

## min_var_single function

### Inputs:
#### base.w: a single base observations weight that is a solution to the constraint
#### Add.Obs.w: the additional observation weight bases that do not affect the expectation
#### ASigAT: the matrix resulting from the quadratic form of A and Sigma
#### method: the method for optimization in the optim() function (default="CG")
#### maxit: passed to the control input of optim() (default=10000)
### Outputs: List of the following:
#### Cvec: the optimal observation weights
#### Variance: the optimal (minimum) variance, based on input Sigma
#### DID.weights: one set of corresponding DID weights

min_var_single <- function(base.w,Add.Obs.w,ASigAT,method="CG", maxit=10000) {
  W_mat <- as.matrix(cbind(as.matrix(base.w, ncol=1),Add.Obs.w))
  colnames(W_mat) <- c("base.w",colnames(Add.Obs.w))
  Inn_mat <- t(W_mat) %*% ASigAT %*% W_mat
  if (ncol(Add.Obs.w)==1) { ## Solved quadratic case for k=1:
    if (Inn_mat[2,2] > 0) {
      Min <- -1*Inn_mat[1,2]/Inn_mat[2,2]
      Var <- Inn_mat[1,1]+2*Min*Inn_mat[1,2]+Min^2*Inn_mat[2,2]
      # Var <- as.numeric(opt_fn(x=Min, Inner=Inn_mat))
    } else {
      stop(simpleError("No Minimum Found."))
    }
  } else {
    Opt <- optim(par=rep(1,ncol(Add.Obs.w)),
                 fn=opt_fn,
                 method=method, control=list(maxit=maxit),
                 Inner=Inn_mat)
    Min <- Opt$par
    Var <- Opt$value
  }
  return(list(Cvec=c(1,Min),
              Variance=Var,
              DID.weights=W_mat %*% matrix(c(1,Min),ncol=1)))
}

# min_var function

### Inputs:
#### SolveOut: The output from the Solve_Assumption function
#### Sigma: The variance matrix to use for variance minimization
#### Drop_ADFT: If TRUE, the ADFT will not be passed through to the output
#### method: the method for optimization in the optim() function (default="CG")
#### maxit: passed to the control input of optim() (default=10000)
### Outputs: List of the following:
#### Cvec: the optimal observation weights
#### Variance: the optimal (minimum) variance
#### DID.weights: one set of corresponding DID weights
#### ADFT: the ADFT object passed through from SolveOut

## A wrapper that outputs the full results of the minimization and can handle a matrix of v's:
## Note the Var that is returned is based on the input Sigma,
### and should generally be used only for relative comparisons of weightings
min_var <- function(SolveOut,
                    Sigma,
                    Drop_ADFT=FALSE,
                    method="CG",maxit=10000) {
  A_mat <- SolveOut$ADFT$A_mat
  solve_obj <- SolveOut$Solve
  if (nrow(Sigma) != ncol(Sigma) | nrow(Sigma) != ncol(A_mat)) {
    stop(simpleError("Sigma must be a square matrix with the same number of columns as A."))
  }

  base.w <- solve_obj$DID.weights %>% dplyr::select(starts_with("w.base")) %>%
    dplyr::rename_all(~stringr::str_replace(.,"^w.base.",""))
  Add.Obs.w <- solve_obj$DID.weights %>% dplyr::select(starts_with("Add.Obs.weights"))
  ASigAT <- A_mat %*% Sigma %*% t(A_mat)

  if (ncol(Add.Obs.w) == 0) {
    print("There is a unique solution to the constraint.")
    res <- list(Cvec=c(1),
                Variance=diag(t(base.w) %*% ASigAT %*% as.matrix(base.w)),
                DID.weights=base.w,
                Obs.weights=solve_obj$Obs.weights %>% dplyr::select(starts_with("ATw.base")))
  } else {
    if (ncol(base.w)==1) { ## A single target estimand
      res1 <- min_var_single(base.w,Add.Obs.w,ASigAT,
                            method=method,maxit=maxit)
      Cvec <- as.matrix(res1$Cvec, ncol=1)
      colnames(Cvec) <- colnames(base.w)
      Var <- as.matrix(res1$Variance, ncol=1)
      colnames(Var) <- colnames(base.w)
      DID.weights <- as.matrix(res1$DID.weights, ncol=1)
      colnames(DID.weights) <- colnames(base.w)
      res <- list(Cvec=Cvec,
                  Variance=Var,
                  DID.weights=DID.weights,
                  Obs.weights=t(A_mat) %*% DID.weights)
    } else { ## Multiple target estimands
      res1 <- apply(base.w, MARGIN=2,
                   FUN=function(col) min_var_single(col,Add.Obs.w,ASigAT,
                                                    method=method,maxit=maxit))
      DID.weights <- do.call(cbind, lapply(res1, FUN=function(x) x$DID.weights))
      colnames(DID.weights) <- colnames(base.w)
      res <- list(Cvec=do.call(cbind, lapply(res1, FUN=function(x) x$Cvec)),
                  Variance=do.call(cbind, lapply(res1, FUN=function(x) x$Variance)),
                  DID.weights=DID.weights,
                  Obs.weights=t(A_mat) %*% DID.weights)
    }
  }
  if (Drop_ADFT) {
    return(res)
  } else {
    return(list(ADFT=SolveOut$ADFT,
                MV=res))
  }
}





#     Min_Res <- apply(v, MARGIN=2,
#                      FUN=function(col) min_var_single()
#     dim_C <- ncol(Add.Obs.w)
#     if (is.vector(v)) {
#       v <- matrix(data=v, ncol=1)
#       ### Make a function out of this as a function of v, use that to vectorize
#       W_mat <- as.matrix(cbind(base.w,Add.Obs.w))
#       Inn_mat <- t(W_mat) %*% A_mat %*% Sigma %*% t(A_mat) %*% W_mat
#       opt_fn <- function(x) {matrix(data=c(1,x), nrow=1) %*% Inn_mat %*% matrix(data=c(1,x), ncol=1)}
#       if (dim_C==1) { ## Solved quadratic case for k=1:
#         if (Inn_mat[2,2] > 0) {
#           Min <- -1*Inn_mat[1,2]/Inn_mat[2,2]
#           Var <- as.numeric(opt_fn(Min))
#         } else {
#           stop(simpleError("No Minimum Found."))
#         }
#       } else {
#         Opt <- optim(par=rep(0,dim_C),
#                      fn=opt_fn,
#                      method="Nelder-Mead")
#         Min <- Opt$par
#         Var <- Opt$value
#       }
#     } else {
#       ## Need to iterate through each v
#     }
#     return(list(Cvec=c(1,Min),
#                 Variance=Var,
#                 DID.weights=W_mat %*% matrix(c(1,Min),ncol=1),
#                 Obs.weights=t(A_mat) %*% W_mat %*% matrix(c(1,Min),ncol=1)))
#   }
# }
