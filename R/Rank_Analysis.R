#######################################
###### File: Rank_Analysis.R ##########
###### Lee Kennedy-Shaffer ############
#######################################

# rank_an function

### Inputs:
#### ADFT_obj: the output from gen_ADFT
#### v: a vector or matrix of estimand weights (e.g., output of create_V)
### Output: List of the following:
#### FT_qr: QR decomposition of F matrix
#### RankAT: Rank of the transpose of A
#### Length_w: Length of the w vector
#### FTv_Ranks: Relevant ranks and W space dimensions for the v's, in order

rank_an <- function(ADFT_obj, v) {
  ## Turn vector v into a matrix if needed:
  if (is.vector(v)) {
    v <- matrix(data=v, ncol=1)
  }

  ## Check dimension compatibility of v and F_mat:
  if (nrow(v) != ncol(ADFT_obj$F_mat)) {
    stop(simpleError("v must be a vector with length equal to the number of columns in F or a matrix of such column vectors."))
  }

  ## Pull key features that don't depend on v:
  RankAT <- (ADFT_obj$N-1)*(ADFT_obj$J-1)
  F_mat <- ADFT_obj$F_mat
  FT_qr <- qr(x=t(F_mat), LAPACK=FALSE)
  Length_w <- (ADFT_obj$N)*(ADFT_obj$N-1)*(ADFT_obj$J)*(ADFT_obj$J-1)/4

  FTv_Rank_Res <- function(v) {
    FTv_Rank <- qr(x=cbind(t(F_mat),v), LAPACK=FALSE)$rank
    if (FTv_Rank > FT_qr$rank) {
      Dim_W <- -1
      Dim_WA <- -1
    } else if (FTv_Rank == FT_qr$rank) {
      Dim_W <- Length_w - FTv_Rank
      Dim_WA <- RankAT - FTv_Rank
      if (Dim_WA < 0) {
        stop(simpleError(paste0("An error has occurred, giving rank(F) = ",FTv_Rank,
                                " > rank(A) = ",RankAT,". Please check inputs.")))
      }
    } else {
      stop(simpleError(paste0("An error has occurred, giving rank(F'|v) = ",FTv_Rank," < rank(F') = ",
                              FT_qr$rank,". Please check inputs.")))
    }
    return(c(FTv_Rank,Dim_W,Dim_WA))
  }

  ## Check ranks across v columns:
  FTv_Ranks <- apply(v, MARGIN=2,
                     FUN=FTv_Rank_Res)
  rownames(FTv_Ranks) <- c("FTv_Rank","Dim_W","Dim_WA")

  if (ncol(FTv_Ranks)==1) {
    if (FTv_Ranks["Dim_W",1] < 0) {
      print(paste0("There are no unique solutions for v"))
    } else if (FTv_Ranks["Dim_W",1]==0) {
      print(paste0("There is one unique solution for v"))
    } else if (FTv_Ranks["Dim_WA",1]==0) {
      print(paste0("There are ",FTv_Ranks["Dim_W",1]," dimensions of DID weight solutions but only one unique observation weight solution for v"))
    } else {
      print(paste0("There are ",FTv_Ranks["Dim_W",1]," dimensions of DID weight solutions and ",FTv_Ranks["Dim_WA",1]," dimensions of observation weight solutions for v"))
    }
  } else {
    for (i in 1:ncol(FTv_Ranks)) {
      if (FTv_Ranks["Dim_W",i] < 0) {
        print(paste0("There are no unique solutions for column ",i))
      } else if (FTv_Ranks["Dim_W",i]==0) {
        print(paste0("There is one unique solution for column ",i))
      } else if (FTv_Ranks["Dim_WA",i]==0) {
        print(paste0("There are ",FTv_Ranks["Dim_W",i]," dimensions of DID weight solutions but only one unique observation weight solution for column ",i))
      } else {
        print(paste0("There are ",FTv_Ranks["Dim_W",i]," dimensions of DID weight solutions and ",FTv_Ranks["Dim_WA",i]," dimensions of observation weight solutions for column ",i))
      }
    }
  }

  return(list(FT_qr=FT_qr,
              RankAT=RankAT,
              Length_w=Length_w,
              FTv_Ranks=FTv_Ranks))
}
