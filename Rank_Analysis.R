#######################################
###### File: Rank_Analysis.R ##########
###### Lee Kennedy-Shaffer ############
###### Created 2024/04/18 #############
###### Updated 2024/04/18 #############
#######################################

rank_an <- function(DFT_obj, v) {
  if (length(v) != dim(DFT_obj$F_mat)[2]) {
    stop(simpleError("v must be a vector with length equal to the number of columns in F"))
  }
  Rank_A <- (DFT_obj$N-1)*(DFT_obj$J-1)
  Rank_F <- qr(DFT_obj$F_mat)$rank
  Rank_Fv <- qr(cbind(t(DFT_obj$F_mat),v))$rank
  Length_w <- (DFT_obj$N)*(DFT_obj$N-1)*(DFT_obj$J)*(DFT_obj$J-1)/4
  if (Rank_Fv > Rank_F) {
    print("The constraint has no solutions.")
    ## Sets dimension variables to negative 1 to indicate no solutions.
    dim_W <- -1
    dim_WA <- -1
  } else if (Rank_Fv == Rank_F) {
    if (Rank_F == Length_w) {
      print("The constraint has a unique solution.")
      ## Sets dimension variables to 0 to indicate unique solution.
      dim_W <- 0
      dim_WA <- 0
    } else if (Rank_F < Length_w) {
      dim_W <- Length_w - Rank_F
      print(paste0("The constraint has an infinite number of solutions, dimension ",dim_W," in terms of the DID estimator weights."))
      dim_WA <- Rank_A - Rank_F
      if (dim_WA == 0) {
        print("There is a unique solution in terms of the observation weights.")
      } else if (dim_WA > 0) {
        print(paste0("The constraint has an infinite number of solutions, dimension ",dim_WA," in terms of the observation weights."))
      } else {
        stop(simpleError(paste0("An error has occurred, giving rank(F) = ",Rank_F," > rank(A) = ",Rank_A,". Please check inputs.")))
      }
    } else {
      stop(simpleError(paste0("An error has occurred, giving rank(F) = ",Rank_F," > dim(w) = ",Length_w,". Please check inputs.")))
    }
  } else {
    stop(simpleError(paste0("An error has occurred, giving rank(F'|v) = ",Rank_Fv," < rank(F') = ",Rank_F,". Please check inputs.")))
  }
  
  return(list(Rank_A=Rank_A,
              Rank_F=Rank_F,
              Rank_Fv=Rank_Fv,
              Length_w=Length_w,
              Dimensions=c(W=dim_W,WA=dim_WA)
              ))
}
