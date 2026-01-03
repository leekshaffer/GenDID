#######################################
###### File: Permutation_Fns.R ########
###### Lee Kennedy-Shaffer ############
#######################################

# Permute_data function

### Inputs:
#### StartTimes: a DF with columns Cluster and StartPd, giving the unique clusters and start periods
#### J: total number of observed periods
#### data: DF with the outcomes, ordered by unit and then period
#### Obs.weights: output of min_var function
### Output: value of the observation weights applied to the data

# Permute_data <- function(StartTimes, J, data, Obs.weights) {
#   ST_P <- tibble(StartPd=rep(sample(StartTimes$StartPd, size=nrow(StartTimes), replace=FALSE), each=J),
#                  Cluster=rep(sample(StartTimes$Cluster, size=nrow(StartTimes), replace=FALSE), each=J))
#   return(t(Obs.weights) %*% as.matrix(cbind(ST_P, data) %>% dplyr::arrange(StartPd,Cluster) %>%
#                                         dplyr::select(-c(Cluster,StartPd))))
# }

# Permute_obs function

## Permute_order helper function

### Inputs:
#### N: total number of units
#### J: total number of observed periods
### Output: randomly permuted unit order

Permute_order <- function(N, J) {
  return((rep(sample.int(N, size=N, replace=FALSE), each=J)-1)*J+rep(1:J, times=N))
}

## Permute_obs main function

### Inputs:
#### Observations: vector of the observed outcomes, ordered by unit then period
#### N: total number of units
#### J: total number of observed periods
#### Obs.weights: the Obs.weights element of the output of min_var function (or NULL)
#### CI.Tx.List: the output of CI_get_Tx_All: a list of ordered, framed treatment effects to test in CIs (or NULL)
#### Drop_Obs: logical; if TRUE, the set of observations will not be included in the output
####  (useful for saving memory when the function will be repeated)
### Output: if Obs.weights is NULL, observations using the randomly permuted unit order
###  Otherwise a list of the following:
#### Obs: observations using the randomly permuted unit order (unless Drop_Obs=TRUE)
#### Ests: the estimate from applying the Obs.weights to the permuted order
#### Additional elements named according to CI.Tx.List if that input is given

Permute_obs <- function(Observations, N, J,
                        Obs.weights=NULL,
                        CI.Tx.List=NULL,
                        Drop_Obs=FALSE) {

  Order <- Permute_order(N, J)
  if (is.null(Obs.weights)) {
    return(Observations[Order,])
  } else {
    Obs <- Observations[Order,]
    Est <- t(Obs.weights) %*% as.matrix(Observations)[Order,,drop=FALSE]
    if (Drop_Obs) {
      out <- list(Ests=Est)
    } else {
      out <- list(Obs=Obs,
                  Ests=Est)
    }
    if (!is.null(CI.Tx.List)) {
      CI.ests <- lapply(X=CI.Tx.List,
             FUN=function(x) Est - t(Obs.weights) %*% as.matrix(x)[Order,,drop=FALSE])
      out <- c(out, CI.ests)
    }
    return(out)
  }
}

# CI helper functions

## CI_get_Tx helper function

### Inputs:
#### Obs_Frame: the full list of clusters and periods, ordered according to Observations
####  (by cluster, earliest switchers first, then by period)
#### Th_All: the Theta$All output from gen_ADFT
#### Th_len: the number of unique treatment effects in the assumption setting (number of rows in Th_All)
#### CI_mat: a vector or matrix of column vectors of treatment effects to test for inclusion in CI, ordered as Th_All
####  If a matrix of multiple columns, they must be the same order as the columns of Observations.
### Output:
#### Data frame of the ordered treatment effects to adjust, corresponding to the full Obs_Frame

CI_get_Tx <- function(Obs_Frame,
                       Th_All, Th_len,
                       CI_mat) {
  ## Converting vector CI_mat to data frame/tibble
    if (is.null(dim(CI_mat))) {
      CI_mat <- tibble(CI_mat)
    }
    ## Fixing length of CI_mat if necessary and giving it theta numbering
    if (nrow(CI_mat) == 1) {
      warning(simpleWarning(paste("CI_mat contains only one value. It will be repeated to", Th_len, "values.")))
      Th_All <- Th_All %>% cbind(CI_mat)
    } else if (nrow(CI_mat) == Th_len) {
      Th_All <- Th_All %>% cbind(CI_mat)
    } else if (nrow(CI_mat) > Th_len) {
      warning(simpleWarning(paste("CI_mat is longer than the number of unique treatment effects. It will be truncated to", Th_len, "values.")))
      Th_All <- Th_All %>% cbind(CI_mat[1:Th_len,,drop=FALSE])
    } else if (nrow(CI_mat) < Th_len) {
      warning(simpleWarning(paste("CI_mat is shorter than the number of unique treatment effects. It will be padded with 0s to", Th_len, "values.")))
      zero_pad <- data.frame(matrix(data=0, nrow=Th_len-nrow(CI_mat), ncol=ncol(CI_mat)))
      names(zero_pad) <- names(CI_mat)
      Th_All <- Th_All %>% cbind(rbind(CI_mat,zero_pad))
    } else {
      stop(simpleError(paste("Error in Dimension/Type of CI_list element",i,". Skipping this element.")))
    }
  ## Add a row for NA theta:
  NA_row <- c(NA,0,0)
  Th_All <- Th_All %>% rbind(NA_row)
  ## Put treatment effect value for each Cl.Num and Pd.Num in Obs_Frame:
  Tx_Ordered <- Obs_Frame %>% left_join(Th_All, by=join_by(Theta)) %>%
    dplyr::select(-c(Cl.Num,Pd.Num,Theta))
  return(Tx_Ordered)
}

## CI_get_Tx_All helper function

### Inputs:
#### ADFT_obj: the gen_ADFT object for the relevant assumption setting
#### CI_list: a (preferably named) list of vector or matrix of column vectors of treatment effects to test for inclusion in CI
### Output:
#### (Named) list of data frames of the ordered treatment effects to adjust, corresponding to the full Obs_Frame

CI_get_Tx_All <- function(ADFT_obj, CI_list) {
  Th_All <- ADFT_obj$Theta$All %>% dplyr::select(Theta)
  Th_len <- nrow(Th_All)
  Obs_Frame <- tibble(Cl.Num=rep(1:N, each=J),
                      Pd.Num=rep(1:J, times=N)) %>%
    left_join(ADFT_obj$Theta$Full %>% dplyr::select(Cl.Num,Pd.Num,Theta),
              by=join_by(Cl.Num,Pd.Num))
  return(lapply(X=CI_list,
                FUN=function(x) CI_get_Tx(Obs_Frame,
                                          Th_All, Th_len,
                                          x)))
}
