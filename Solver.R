#######################################
###### File: Solver.R #################
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/18 #############
###### Updated 2024/04/18 #############
#######################################

## Note v can be a vector to solve for one target,
### or it can be a matrix of column vectors to solve for.
solve_WA <- function(DFT_obj,A_mat,v,DID_full=FALSE) {
  ## If v is a vector, put it into matrix form:
  if (is.vector(v)) {
    v <- matrix(data=v, ncol=1)
  }
  
  ## Get the key info from A and F
  F_mat <- DFT_obj$F_mat
  F_qr <- qr(x=t(F_mat))
  RankAT <- (DFT_obj$N-1)*(DFT_obj$J-1)
  
  ## Check Rank Conditions for each v
  FTv.Check <- apply(v, MARGIN=2,
        FUN=function(col) qr(x=cbind(t(F_mat),col))$rank > F_qr$rank)
  if (sum(FTv.Check)==0) {
  } else if (sum(FTv.Check)==ncol(v)) {
    stop(simpleError("No columns of v have solutions."))
  } else {
    warning(simpleWarning(paste0("The following columns of v have no solutions and were dropped: ",
                                 paste((1:ncol(v))[FTv.Check], collapse=", "),
                                 ". The remaining columns have been shifted accordingly. Interpret results accordingly.")))
    v <- v[,!FTv.Check,drop=FALSE]
  }
  
  ## Check for column of all zero's in F and corresponding 0 row in v
  ZeroCols <- apply(F_mat, MARGIN=2,
                    FUN=function(col) sum(col==0) == nrow(F_mat))
  if (sum(ZeroCols) > 0) {
    if (sum(v[ZeroCols,]==0) != ncol(v)) {
      stop(simpleError("v has values that cannot be achieved with this F matrix."))
    } else {
      F_mat <- F_mat[,!ZeroCols,drop=FALSE]
      F_qr <- qr(x=t(F_mat))
      v <- v[!ZeroCols,,drop=FALSE]
    }
  }

  
  if (nrow(v) != ncol(F_mat)) {
    stop(simpleError("v must be a vector with length corresponding to the number of columns in F_mat, or a matrix of such column vectors."))
  } else {
    ## find a single solution for w:
    w <- qr.solve(a=F_qr, b=v)
    DID.weights <- data.frame(w.base=w)
    
    ## find a single solution for weights for observations via A' * w:
    Obs.weights <- data.frame(ATw.base=t(A_mat) %*% w)
    
    if (F_qr$rank==nrow(w)) { ## Prints note for single solution and jumps to return
      print("There is a unique solution for the DID estimator weights.")
    } else { ## If there are non-unique DID estimator weights
      ## Use SVD of A' to find basis for its kernel
      AT_svd <- svd(x=t(A_mat),
                    nu=0,
                    nv=nrow(A_mat))
      kerAT_basis <- AT_svd$v[,(RankAT+1):ncol(AT_svd$v),drop=FALSE]
      ### Re-normalize the basis to sum (in abs. value) to 1 and smooth near-zeros:
      kerAT_norm <- apply(kerAT_basis, MARGIN=2,
                          FUN=function(col) ifelse(abs(col) < .Machine$double.eps, 0, col/sum(abs(col))))
      
      ## If rank(F) < rank(A), get basis of ker(F') that is orthogonal to ker(A')
      if (F_qr$rank < RankAT) {
        ## Use SVD of F' to find basis for its kernel
        FT_svd <- svd(x=t(F_mat),
                      nu=0, nv=nrow(F_mat))
        kerFT_basis <- FT_svd$v[,(F_qr$rank+1):ncol(FT_svd$v), drop=FALSE]
        
        ## Use QR decomposition on (A' basis | F' basis) to get Q, 
        ### whose first nullity(F') columns are an orthonormal basis for F',
        ### the first nullity(A') columns of which are an orthonormal basis for A',
        ### so the nullity(F')-nullity(A') columns between are an orthonormal basis for F'\A'
        kerFT_only <- qr.Q(qr(cbind(kerAT_basis,kerFT_basis)))[,(dim(kerAT_basis)[2]+1):(dim(kerFT_basis)[2]), drop=FALSE]
        ### Re-normalize them to sum (in abs. value) to 1 and smooth near-zeros:
        kerFT_norm <- apply(kerFT_only, MARGIN=2,
                            FUN=function(col) ifelse(abs(col) < .Machine$double.eps, 0, col/sum(abs(col))))
        
        ## Get observation weights of these F'\A' basis vectors
        ATw.weights <- t(A_mat) %*% kerFT_only
        ### Re-normalize them to sum (in abs. value) to 1 and smooth near-zeros:
        # ATw_norm <- apply(ATw.weights, MARGIN=2,
        #                     FUN=function(col) ifelse(abs(col) < .Machine$double.eps, 0, col/sum(abs(col))))
        ## Use normalized DID weights to get observation weights:
        ATw_norm <- apply(t(A_mat) %*% kerFT_norm, MARGIN=2,
                          FUN=function(col) ifelse(abs(col) < .Machine$double.eps, 0, col))
        
        ## Append F'\A' basis vector weights to Obs.weights and DID.weights as "Add.Obs.weights"
        Obs.weights <- cbind(Obs.weights,Add.Obs.weights=ATw_norm)
        DID.weights <- cbind(DID.weights,Add.Obs.weights=kerFT_norm)
      } else if (F_qr$rank > RankAT) {
        stop(simpleError("F has greater rank than A. Please check inputs for accuracy."))
      }
      
      if (DID_full) { ## Add DID-estimator-only weights (that do not affect observation weights) if requested
        print("Note: the weights labelled Add.DID.weights affect the estimator weights but not the observation weights.")
        DID.weights <- cbind(DID.weights,Add.DID.weights=kerAT_norm)
      } else { ## Otherwise print a note
        print("Note: returning only a single solution for the DID estimator weights, althoughs others may exist that yield equivalent observation weights.")
      }
    }
    return(list(DID.weights=DID.weights,
                Obs.weights=Obs.weights))
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
