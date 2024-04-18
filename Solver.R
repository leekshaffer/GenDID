#######################################
###### File: Solver.R #################
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/18 #############
###### Updated 2024/04/18 #############
#######################################

## Method for getting the solutions when dim(WA) = 0 and dim(W) >= 0
## Should use only when Rank(A) = Rank(F) = Rank(F'|v):
## Note v can be a vector to solve for one target,
### or it can be a matrix of column vectors to solve for.
solve_WA0 <- function(DFT_obj,A_mat,v,DID_full=TRUE) {
  if (is.vector(v)) {
    v <- matrix(data=v, ncol=1)
  }
  
  if (dim(v)[1] != dim(DFT_obj$F_mat)[2]) {
    stop(simpleError("v must be a vector with length corresponding to the number of columns in F_mat, or a matrix of such column vectors."))
  } else {
    ## find a single solution for w:
    F_qr <- qr(x=t(DFT_obj$F_mat))
    w <- qr.solve(a=F_qr, b=v)
    
    ## find unique weights for observations:
    ATw <- t(A_mat) %*% w
    
    ## find full DID weights solution set if requested:
    if (DID_full) {
      if (F_qr$rank==dim(w)[1]) {
        print("There is a unique solution for the DID estimator weights.")
        D_aug <- cbind(Weights=w,DFT_obj$D_aug)
      } else {
        AT_svd <- svd(x=t(A_mat),
                      nu=0,
                      nv=dim(A_mat)[1])
        RankAT <- (DFT_obj$N-1)*(DFT_obj$J-1)
        kerAT_basis <- AT_svd$v[,(RankAT+1):ncol(AT_svd$v),drop=FALSE]
        D_aug <- cbind(Weights=w, Addl.Weight=kerAT_basis,
                       DFT_obj$D_aug)
      }
    } else {
      print("Returning only a single solution for the DID estimator weights.")
      D_aug <- cbind(Weights=w,DFT_obj$D_aug)
    }
    return(list(DID.weights=w,
                Obs.weights=ATw,
                D_aug=D_aug))
  }
}

solve_WA <- function(DFT_obj,A_mat,v,DID_full=TRUE) {
  if (is.vector(v)) {
    v <- matrix(data=v, ncol=1)
  }
  
  if (dim(v)[1] != dim(DFT_obj$F_mat)[2]) {
    stop(simpleError("v must be a vector with length corresponding to the number of columns in F_mat, or a matrix of such column vectors."))
  } else {
    ## find a single solution for w:
    F_qr <- qr(x=t(DFT_obj$F_mat))
    w <- qr.solve(a=F_qr, b=v)
    DID.weights <- data.frame(w.base=w)
    
    ## find a single solution for weights for observations:
    Obs.weights <- data.frame(ATw.base=t(A_mat) %*% w)
    
    if (F_qr$rank==dim(w)[1]) {
      print("There is a unique solution for the DID estimator weights.")
    } else {
      AT_svd <- svd(x=t(A_mat),
                    nu=0,
                    nv=dim(A_mat)[1])
      RankAT <- (DFT_obj$N-1)*(DFT_obj$J-1)
      kerAT_basis <- AT_svd$v[,(RankAT+1):ncol(AT_svd$v),drop=FALSE]
      
      if (F_qr$rank == RankAT) {
        if (DID_full) {
          DID.weights <- cbind(DID.weights,Add.DID.weights=kerAT_basis)
        } else {
          print("Note: returning only a single solution for the DID estimator weights, althoughs others may exist that yield equivalent observation weights.")
        }
      } else if (F_qr$rank < RankAT) {
        FT_svd <- svd(x=t(DFT_obj$F_mat),
                      nu=0, nv=dim(DFT_obj$F_mat)[1])
        kerFT_basis <- FT_svd$v[,(F_qr$rank+1):ncol(FT_svd$v), drop=FALSE]
        kerFT_only <- qr.Q(qr(cbind(kerAT_basis,kerFT_basis)))[,(dim(kerAT_basis)[2]+1):(dim(kerFT_basis)[2]), drop=FALSE]
        ATw.weights <- t(A_mat) %*% kerFT_only
        Obs.weights <- cbind(Obs.weights,Add.Obs.weights=ATw.weights)
        DID.weights <- cbind(DID.weights,Add.Obs.weights=kerFT_only)
        
        if (DID_full) {
         DID.weights <- cbind(DID.weights,Add.DID.weights=kerAT_basis)
        } else {
          print("Note: returning only a single solution for the DID estimator weights, althoughs others may exist that yield equivalent observation weights.")
        }
      } else {
        print("Error: Rank(F) calculated to be greater than Rank(A). Please check inputs.")
      }
    }
    return(list(DID.weights=DID.weights,
                Obs.weights=Obs.weights,
                D_aug=D_aug))
  }
}

# solve_svd <- function(DFT_obj,A_mat,v) {
#   FTsvd <- svd(x=t(DFT_obj$F_mat),
#                nu=dim(DFT_obj$F_mat)[2],
#                nv=dim(DFT_obj$F_mat)[1])
#   ATsvd <- svd(x=t(A_mat),
#                nu=dim(A_mat)[2],
#                nv=dim(A_mat)[1])
# }
# 
# require(MASS)
# 
# solve_gen <- function(DFT_obj,A_mat,v) {
#   MPg <- MASS::ginv(t(DFT_obj$F_mat))
#   Soln1 <- MPg %*% v
#   Free <- diag(x=1,nrow=dim(MPg)[1]) - MPg %*% t(DFT_obj$F_mat)
# }
