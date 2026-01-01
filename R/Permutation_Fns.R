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
#### Obs.weights: output of min_var function (or NULL)
### Output: if Obs.weights is NULL, observations using the randomly permuted unit order
###  Otherwise a list of the following:
#### Obs: observations using the randomly permuted unit order
#### Ests: the estimate from applying the Obs.weights to the permuted order

Permute_obs <- function(Observations, N, J, Obs.weights=NULL) {
  Order <- Permute_order(N, J)
  if (is.null(Obs.weights)) {
    return(Observations[Order,])
  } else {
    return(list(Obs=Observations[Order,],
                Ests=t(Obs.weights) %*% as.matrix(Observations)[Order,,drop=FALSE]))
  }
}
