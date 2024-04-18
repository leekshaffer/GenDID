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
  } else if (Rank_Fv == Rank_F) {
    if (Rank_F == Length_w) {
      print("The constraint has a unique solution.")
    } else if (Rank_F < Length_w) {
      print(paste0("The constraint has an infinite number of solutions, dimension ",Length_w-Rank_F," in terms of the DID estimator weights."))
      diff_A_F <- Rank_A - Rank_F
      if (diff_A_F == 0) {
        print("There is a unique solution in terms of the observation weights.")
      } else if (diff_A_F > 0) {
        print(paste0("The constraint has an infinite number of solutions, dimension ",diff_A_F," in terms of the observation weights."))
      } else {
        print("An error has occurred, giving rank(F) > rank(A). Please check inputs.")
      }
    } else {
      print("An error has occurred, giving rank(F) > dim(w). Please check inputs.")
    }
  } else {
    print("An error has occurred, giving rank(F'|v) < rank(F'). Please check inputs.")
  }
  
  return(list(Rank_A=Rank_A,
              Rank_F=Rank_F,
              Rank_Fv=Rank_Fv,
              Length_w=Length_w))
}
