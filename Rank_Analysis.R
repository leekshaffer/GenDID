#######################################
###### File: Rank_Analysis.R ##########
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/18 #############
###### Updated 2024/04/18 #############
#######################################

rank_an <- function(DFT_obj, v) {
  ## Turn vector v into a matrix if needed:
  if (is.vector(v)) {
    v <- matrix(data=v, ncol=1)
  }
  
  ## Check dimension compatibility of v and F_mat:
  if (nrow(v) != ncol(DFT_obj$F_mat)) {
    stop(simpleError("v must be a vector with length equal to the number of columns in F or a matrix of such column vectors."))
  }
  
  ## Pull key features that don't depend on v:
  RankAT <- (DFT_obj$N-1)*(DFT_obj$J-1)
  F_mat <- DFT_obj$F_mat
  FT_qr <- qr(x=t(F_mat))
  Length_w <- (DFT_obj$N)*(DFT_obj$N-1)*(DFT_obj$J)*(DFT_obj$J-1)/4
  
  FTv_Rank_Res <- function(v) {
    FTv_Rank <- qr(x=cbind(t(F_mat),col))$rank
    if (FTv_Rank > FT_qr$rank) {
      Dim_W <- -1
      Dim_WA <- -1
    } else if (FTv_Rank == FT_qr$rank) {
      Dim_W <- FTv_Rank - Length_w
      Dim_WA <- RankAT - FTv_Rank
      if (Dim_WA < 0) {
        stop(simpleError(paste0("An error has occurred, giving rank(F) = ",FTv_Rank,
                                " > rank(A) = ",RankAT,". Please check inputs.")))
      }
    } else {
      stop(simpleError(paste0("An error has occurred, giving rank(F'|v) = ",FTv_Rank," < rank(F') = ",
                              FT_qr$rank,". Please check inputs.")))
    }
    return(rbind(FTv_Rank=FTv_Rank,
                 Dim_W=Dim_W,
                 Dim_WA=Dim_WA))
  }

  ## Check ranks across v columns:
  FTv_Ranks <- apply(v, MARGIN=2,
                     FUN=##apply above function
  )
  
  
  # if (Rank_Fv > Rank_F) {
  #   print("The constraint has no solutions.")
  #   ## Sets dimension variables to negative 1 to indicate no solutions.
  #   dim_W <- -1
  #   dim_WA <- -1
  # } else if (Rank_Fv == Rank_F) {
  #   if (Rank_F == Length_w) {
  #     print("The constraint has a unique solution.")
  #     ## Sets dimension variables to 0 to indicate unique solution.
  #     dim_W <- 0
  #     dim_WA <- 0
  #   } else if (Rank_F < Length_w) {
  #     dim_W <- Length_w - Rank_F
  #     print(paste0("The constraint has an infinite number of solutions, dimension ",dim_W," in terms of the DID estimator weights."))
  #     dim_WA <- Rank_A - Rank_F
  #     if (dim_WA == 0) {
  #       print("There is a unique solution in terms of the observation weights.")
  #     } else if (dim_WA > 0) {
  #       print(paste0("The constraint has an infinite number of solutions, dimension ",dim_WA," in terms of the observation weights."))
  #     } else {
  #       
  #     }
  #   } else {
  #     stop(simpleError(paste0("An error has occurred, giving rank(F) = ",Rank_F," > dim(w) = ",Length_w,". Please check inputs.")))
  #   }
  # } else {
  #   
  # }
  
  ## Update return to keep the ranks and stuff.
  return(list(Rank_A=Rank_A,
              Rank_F=Rank_F,
              Rank_Fv=Rank_Fv,
              Length_w=Length_w,
              Dimensions=c(W=dim_W,WA=dim_WA)
              ))
}
