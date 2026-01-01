#######################################
###### File: Solver.R #################
###### Lee Kennedy-Shaffer ############
#######################################

# solve_WA function

### Inputs:
#### ADFT_obj: the output from gen_ADFT
#### v: a vector or matrix of estimand weights (e.g., output of create_V)
#### rank_obj: the output from rank_an (if available)
#### DID_full: Boolean:
####  TRUE adds weight vectors that only change DID decomposition, not observation weights;
####  FALSE omits these
### Output: List of the following:
#### DID.weights: basis vectors for the DID estimators weights that solve the constraints
#### Obs.weights: basis vectors for the observation weights that solve the constraints

solve_WA <- function(
    ADFT_obj,
    v,
    rank_obj = NULL,
    DID_full = TRUE
) {
  ## If v is a vector, put it into matrix form:
  if (is.vector(v)) {
    v <- matrix(data = v, ncol = 1)
  }

  F_mat <- ADFT_obj$F_mat

  A_mat <- ADFT_obj$A_mat

  if (!is.null(rank_obj)) { # if given a rank object, use it:
    if (ncol(v) != ncol(rank_obj$FTv_Ranks)) {
      stop("v and rank_obj do not come from the same dimension of v.")
    }

    FT_qr <- rank_obj$FT_qr
    FT_rank <- FT_qr$rank
    RankAT <- rank_obj$RankAT

    # Check for solvability
    no_solution <- rank_obj$FTv_Ranks["Dim_W", ] < 0

    if (all(no_solution)) {
      stop("No columns of v have solutions.")
    } else if (any(no_solution)) {
      warning(paste0(
        "The following columns of v have no solutions and were dropped: ",
        paste(which(no_solution), collapse = ", "),
        ". The remaining columns have been shifted accordingly. Interpret results accordingly."
      ))
      v <- v[, !no_solution, drop = FALSE]
    }
  } else { # otherwise (rank_obj is missing)
    ## Get the key info from A and F
    RankAT <- (ADFT_obj$N - 1) * (ADFT_obj$J - 1)
    FT_qr <- qr(x = t(F_mat), LAPACK = FALSE)
    FT_rank <- FT_qr$rank

    ## Check Rank Conditions for each v
    no_solution <- apply(v, 2,
                      function(col) qr(cbind(t(F_mat), col))$rank > FT_rank)
    if (all(no_solution)) {
      stop("No columns of v have solutions.")
    } else if (any(no_solution)) {
      warning(paste0(
        "The following columns of v have no solutions and were dropped: ",
        paste((1:ncol(v))[FTv.Check], collapse = ", "),
        ". The remaining columns have been shifted accordingly. Interpret results accordingly."
      ))
      v <- v[, !FTv.Check, drop = FALSE]
    }
  }

  if (FT_rank < ncol(F_mat)) { ## Rank-deficient F:

    # Identify columns in F_mat that are entirely zero
    ZeroCols <- colSums(F_mat == 0) == nrow(F_mat)

    if (any(ZeroCols)) {
      if (sum(v[ZeroCols, ] == 0) != ncol(v)) {
        stop("v has values that cannot be achieved with this F matrix.")
      }
      F_mat <- F_mat[, !ZeroCols, drop = FALSE]
      FT_qr <- qr(x = t(F_mat), LAPACK = FALSE)

      v <- v[!ZeroCols, , drop = FALSE]
    }

    # If still rank-deficient, find linear dependent column in F:
    if (FT_qr$rank < ncol(F_mat)) {
      RankCols <- sapply(
        1:ncol(F_mat),
        function(i) qr(F_mat[, 1:i])$rank
      )
      RankPrev <- c(RankCols[1], RankCols[2:length(RankCols)] - RankCols[1:(length(RankCols) - 1)])
      KeepCols <- (1:length(RankPrev))[RankPrev == 1]
      F_mat <- F_mat[, KeepCols, drop = FALSE]
      FT_qr <- qr(t(F_mat), LAPACK = FALSE)
      v <- v[KeepCols, , drop = FALSE]
    }
  }

  if (nrow(v) != ncol(F_mat)) {
    stop("v must be a vector with length corresponding to the number of columns in F_mat, or a matrix of such column vectors.")
  }
  ## find a single solution for w:
  w <- qr.solve(FT_qr, v)
  if (is.null(colnames(w))) {
    colnames(w) <- 1:ncol(w)
  }

  DID.weights <- data.frame(w)
  colnames(DID.weights) <- paste("w.base", colnames(w), sep = ".")

  ## find a single solution for weights for observations via A' * w:
  Obs.weights <- data.frame(t(A_mat) %*% w)
  colnames(Obs.weights) <- paste("ATw.base", colnames(w), sep = ".")

  if (FT_rank == nrow(w)) { ## Prints note for single solution and jumps to return
    print("There is a unique solution for the DID estimator weights.")
    return(
      list(
        DID.weights = DID.weights,
        Obs.weights = Obs.weights
      ))
  } ## If there are non-unique DID estimator weights
  ## Use SVD of A' to find basis for its kernel
  AT_svd <- svd(
    t(A_mat),
    nu = 0,
    nv = nrow(A_mat)
  )
  kerAT_basis <- AT_svd$v[, (RankAT + 1):ncol(AT_svd$v), drop = FALSE]
  ### Re-normalize the basis to sum (in abs. value) to 1 and smooth near-zeros:
  kerAT_norm <- apply(kerAT_basis,
                      MARGIN = 2,
                      FUN = function(col) ifelse(abs(col) < .Machine$double.eps, 0, col / sum(abs(col)))
  )

  ## If rank(F) < rank(A), get basis of ker(F') that is orthogonal to ker(A')
  if (FT_rank < RankAT) {
    ## Use SVD of F' to find basis for its kernel
    FT_svd <- svd(
      t(F_mat),
      nu = 0,
      nv = nrow(F_mat)
    )
    kerFT_basis <- FT_svd$v[, (FT_rank + 1):ncol(FT_svd$v), drop = FALSE]

    ## Use QR decomposition on (A' basis | F' basis) to get Q,
    ### whose first nullity(F') columns are an orthonormal basis for F',
    ### the first nullity(A') columns of which are an orthonormal basis for A',
    ### so the nullity(F')-nullity(A') columns between are an orthonormal basis for F'\A'
    kerFT_only <- qr.Q(qr(cbind(kerAT_basis, kerFT_basis), LAPACK = FALSE))[, (dim(kerAT_basis)[2] + 1):(dim(kerFT_basis)[2]), drop = FALSE]
    ### Re-normalize them to sum (in abs. value) to 1 and smooth near-zeros:
    kerFT_norm <- apply(kerFT_only,
                        2,
                        function(col) ifelse(abs(col) < .Machine$double.eps, 0, col / sum(abs(col)))
    )

    ## Get observation weights of these F'\A' basis vectors
    ATw.weights <- t(A_mat) %*% kerFT_only
    ### Re-normalize them to sum (in abs. value) to 1 and smooth near-zeros:
    # ATw_norm <- apply(ATw.weights, MARGIN=2,
    #                     FUN=function(col) ifelse(abs(col) < .Machine$double.eps, 0, col/sum(abs(col))))
    ## Use normalized DID weights to get observation weights:
    ATw_norm <- apply(t(A_mat) %*% kerFT_norm,
                      2,
                      function(col) ifelse(abs(col) < .Machine$double.eps, 0, col)
    )

    ## Append F'\A' basis vector weights to Obs.weights and DID.weights as "Add.Obs.weights"
    Obs.weights <- cbind(Obs.weights, Add.Obs.weights = ATw_norm)
    DID.weights <- cbind(DID.weights, Add.Obs.weights = kerFT_norm)
  } else if (FT_rank > RankAT) {
    stop("F has greater rank than A. Please check inputs for accuracy.")
  }

  if (DID_full) { ## Add DID-estimator-only weights (that do not affect observation weights) if requested
    print("Note: the weights labelled Add.DID.weights affect the estimator weights but not the observation weights.")
    DID.weights <- cbind(DID.weights, Add.DID.weights = kerAT_norm)
  } else { ## Otherwise print a note
    print("Note: returning only a single solution for the DID estimator weights, althoughs others may exist that yield equivalent observation weights.")
  }

  return(
    list(
      DID.weights = DID.weights,
      Obs.weights = Obs.weights
    ))
}
