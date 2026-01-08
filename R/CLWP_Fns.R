#######################################
###### File: CLWP_Fns.R ###############
#######################################

## Functions to conduct a CLWP Analysis

## See Voldal et al., Statist Med, 2024, https://doi.org/10.1002/sim.10120 ##
## and its supplemental material for details on these methods and the ##
## underlying code. If using these methods, please refer to and cite ##
## Voldal et al. 2024 ##

## Load packages
require(bbmle)

## Helper function: loglik_vertical_CML

### Inputs:
#### theta: a value for the treatment effect
#### sigma: a value for the standard deviation of the vertical contrasts
#### my.data.full: a dataset in the format produced by swCRTdesign
### Output: negative log likelihood

loglik_vertical_CML <- function(theta, sigma, my.data.full){
  #First, we need to go from a data frame of cluster-period observations to a data frame with
  #every vertical contrast between treatment and control cluster-periods.
  my.data.combos <- data.frame(treated=double(),control=double())
  for(j in unique(my.data.full$Period)){
    treated_list.j <- my.data.full$Summ[my.data.full$Period == j & my.data.full$Interv == 1]
    control_list.j <- my.data.full$Summ[my.data.full$Period == j & my.data.full$Interv == 0]
    combos.j <- expand.grid(treated=treated_list.j,control=control_list.j)
    my.data.combos <- rbind(my.data.combos,combos.j)
  }
  my.data.combos$difference <- my.data.combos$treated-my.data.combos$control
  R=suppressWarnings(dnorm(x=my.data.combos$difference,mean=theta,sd=sigma))
  return(-sum(log(R)))
}

## Helper function: get_u_per_contrast

### Inputs:
#### fitted_params: fitted parameters from the vertical CL, e.g. my.fitted.cl@coef, a vector of (theta, sigma)
#### difference: one vertical difference
### Output: u (composite score function) for one vertical difference

get_u_per_contrast <- function(fitted_params,difference){
  theta_cl <- fitted_params["theta"]
  sigmasq_cl <- fitted_params["sigma"]^2#Note: going from sd to var

  u_1 <- (1/sigmasq_cl)*(difference-theta_cl)

  u_2 <- -(1/2)*(1/sigmasq_cl)+(1/2)*(1/sigmasq_cl^2)*(difference-theta_cl)^2

  return(c(u_1,u_2))

}

## Helper function: get_triangledown_u_per_contrast

### Input:
#### fitted_params: fitted parameters from the vertical CL, e.g. my.fitted.cl@coef, a vector of (theta, sigma)
#### difference: one vertical difference
### Output: a matrix of the second derivatives (derivative of the score) for one vertical difference

get_triangledown_u_per_contrast <- function(fitted_params,difference){
  theta_cl <- fitted_params["theta"]
  sigmasq_cl <- fitted_params["sigma"]^2#Note: going from sd to var

  block_1_1 <- -1/sigmasq_cl

  block_1_2 <- -1/sigmasq_cl^2*(difference-theta_cl)

  block_2_2 <- (1/2)*(1/sigmasq_cl^2)-(1/sigmasq_cl^3)*(difference-theta_cl)^2

  row_1 <- c(block_1_1,block_1_2)
  row_2 <- c(block_1_2,block_2_2)

  full_matrix <- rbind(row_1,row_2)

  return(full_matrix)

}

## Helper function: loglik_vertical_CML_adj

### Inputs:
#### theta: a value for the treatment effect
#### sigma: a value for the standard deviation of the vertical contrasts, after adjusting for baseline
#### betabase: a value for the coefficient for the baseline difference
#### my.data.full: a dataset in the format produced by swCRTdesign
### Output: negative log likelihood

loglik_vertical_CML_adj <- function(theta, sigma, betabase, my.data.full){
  #First, we need to go from a data frame of cluster-period observations to a data frame with
  #every vertical contrast between treatment and control cluster-periods.
  my.data.combos <- data.frame(treated=double(),control=double())
  my.data.base.combos <- data.frame(treated=double(),control=double())
  for(j in unique(my.data.full$Period)){
    treated_list.j <- my.data.full$Summ[my.data.full$Period == j & my.data.full$Interv == 1]
    control_list.j <- my.data.full$Summ[my.data.full$Period == j & my.data.full$Interv == 0]
    combos.j <- expand.grid(treated=treated_list.j,control=control_list.j)
    my.data.combos <- rbind(my.data.combos,combos.j)
    #Now do the same thing in the same order, but get the baseline measurements
    #Note: assuming here that the first period is all-control
    #(so won't be contributing to vertical contrasts, and don't need to account for differing trt status)
    treated_list_base.j <- my.data.full$Summ[my.data.full$Period == 1 &
                                                   my.data.full$Cluster %in% my.data.full$Cluster[my.data.full$Period == j & my.data.full$Interv == 1]]
    control_list_base.j <- my.data.full$Summ[my.data.full$Period == 1 &
                                                   my.data.full$Cluster %in% my.data.full$Cluster[my.data.full$Period == j & my.data.full$Interv == 0]]
    combos.base.j <- expand.grid(treated=treated_list_base.j,control=control_list_base.j)
    my.data.base.combos <- rbind(my.data.base.combos,combos.base.j)
  }
  my.data.combos$difference <- my.data.combos$treated-my.data.combos$control
  my.data.base.combos$difference <- my.data.base.combos$treated-my.data.base.combos$control
  R=suppressWarnings(dnorm(x=my.data.combos$difference,mean=theta+betabase*my.data.base.combos$difference,sd=sigma))
  return(-sum(log(R)))
}

## Helper function: get_u_per_contrast_adj

### Inputs:
#### fitted_params: fitted parameters from the vertical CL, e.g. my.fitted.cl@coef, a vector of (theta, sigma, betabase)
#### difference: one vertical difference
#### difference_base: the analogous difference at baseline
### Output: u (composite score function) for one vertical difference

get_u_per_contrast_adj <- function(fitted_params,difference,difference_base){
  theta_cl <- fitted_params["theta"]
  sigmasq_cl <- fitted_params["sigma"]^2#Note: going from sd to var
  betabase_cl <- fitted_params["betabase"]

  u_1 <- (1/sigmasq_cl)*(difference-(theta_cl+betabase_cl*difference_base))

  u_2 <- -(1/2)*(1/sigmasq_cl)+(1/2)*(1/sigmasq_cl^2)*(difference-(theta_cl+betabase_cl*difference_base))^2

  u_3 <- (1/sigmasq_cl)*(difference-(theta_cl+betabase_cl*difference_base))*difference_base


  return(c(u_1,u_2,u_3))

}

## Helper function: get_triangledown_u_per_contrast_adj

### Inputs:
#### fitted_params: fitted parameters from the vertical CL, e.g. my.fitted.cl@coef, a vector of (theta, sigma, betabase)
#### difference: one vertical difference
#### difference_base: the analogous difference at baseline
### Output: a matrix of the second derivatives (derivative of the score) for one vertical difference

get_triangledown_u_per_contrast_adj <- function(fitted_params,difference,difference_base){
  theta_cl <- fitted_params["theta"]
  sigmasq_cl <- fitted_params["sigma"]^2#Note: going from sd to var
  betabase_cl <- fitted_params["betabase"]

  block_1_1 <- -1/sigmasq_cl

  block_1_2 <- -1/sigmasq_cl^2*(difference-(theta_cl+betabase_cl*difference_base))

  block_1_3 <- -difference_base/sigmasq_cl


  block_2_2 <- (1/2)*(1/sigmasq_cl^2)-(1/sigmasq_cl^3)*(difference-(theta_cl+betabase_cl*difference_base))^2

  block_2_3 <- -difference_base/sigmasq_cl^2*(difference-(theta_cl+betabase_cl*difference_base))

  block_3_3 <- -difference_base^2/sigmasq_cl

  row_1 <- c(block_1_1,block_1_2,block_1_3)
  row_2 <- c(block_1_2,block_2_2,block_2_3)
  row_3 <- c(block_1_3,block_2_3,block_3_3)

  full_matrix <- rbind(row_1,row_2,row_3)

  return(full_matrix)

}

## Main CLWP_fit function

### using starting values from GEE
### Note - I did do some experimenting, and it doesn't seem to be very sensitive to starting values.

### Inputs:
#### my.data: data to analyze
#### start_theta: starting value for theta, e.g., summary(gee.independence)$coef["tx.var","Estimate"]
#### start_sigma: starting value for sigma, e.g., sqrt(2*var(gee.independence$residuals))
#### N:
### Output: List of the following:
#### CLWP_est: estimated effect
#### CLWP_se: SE of estimate
#### CLWP_pval: p-value of estimate

CLWP_fit <- function(my.data, start_theta, start_sigma, N) {
  cml.vertical <-mle2(function(theta,sigma){loglik_vertical_CML(theta=theta, sigma=sigma,my.data.full=my.data)},
                      start=list(theta=start_theta, sigma=start_sigma))

  #Recreate data frame of all the differences
  my.data.combos <- data.frame(treated=double(),control=double(),cluster_trt=double(),cluster_ctrl=double())
  for(j in unique(my.data$Period)){
    treated_list.j <- my.data$Summ[my.data$Period == j & my.data$Interv == 1]
    control_list.j <- my.data$Summ[my.data$Period == j & my.data$Interv == 0]
    combos.j <- expand.grid(treated=treated_list.j,control=control_list.j)
    #get cluster ID's for each contrast to pull later
    cluster_treated_list.j <- my.data$Cluster[my.data$Period == j & my.data$Interv == 1]
    cluster_control_list.j <- my.data$Cluster[my.data$Period == j & my.data$Interv == 0]
    cluster_combos.j <- expand.grid(cluster_trt=cluster_treated_list.j,cluster_ctrl=cluster_control_list.j)
    my.data.combos <- rbind(my.data.combos,cbind(combos.j,cluster_combos.j))
  }
  my.data.combos$difference <- my.data.combos$treated-my.data.combos$control


  #Set up empty matrices to sum
  triangledown_u_sum <- matrix(0,nrow=2,ncol=2)#2x2 for the two parameters
  u_sq_sum <- matrix(0,nrow=2,ncol=2)


  for(i in 1:nrow(my.data.combos)){
    triangledown_u_sum <- triangledown_u_sum +get_triangledown_u_per_contrast(fitted_params=cml.vertical@coef,
                                                                              difference=my.data.combos$difference[i])
  }

  #To do J, need to sum u's within clusters, then square, THEN sum over clusters
  for(i in unique(my.data$Cluster)){
    u.i <- c(0,0)#initializing score vector (1x2)
    contrast_involves_i <- my.data.combos$cluster_trt == i | my.data.combos$cluster_ctrl == i
    my.data.combos.i <- my.data.combos[contrast_involves_i,]
    for(k in 1:nrow(my.data.combos.i)){#loop through all the contrasts that involve this cluster
      u.ik <- get_u_per_contrast(fitted_params=cml.vertical@coef,difference=my.data.combos.i$difference[k])
      u.i <- u.i+u.ik
    }
    u_sq_sum <- u_sq_sum + u.i %*% t(u.i)
  }


  my_H <- - triangledown_u_sum/N
  my_J <- u_sq_sum/N
  my_G <- my_H %*% solve(my_J) %*% my_H
  my_variance_matrix <- solve(my_G)

  return(c(CLWP_est=cml.vertical@coef["theta"],
           CLWP_se=((1/N)*my_variance_matrix[1,1])^.5,
           CLWP_pval=2 * pnorm(abs((cml.vertical@coef["theta"])/(((1/N)*my_variance_matrix[1,1])^.5)),
                               lower.tail = FALSE)))
}

## Main CLWPA_fit function

### using same starting values from GEE as the unadjusted version
### the starting standard deviation will be a little too large for the adjusted version
### based on informal experimenting, starting vals too big are fine (too small are riskier).
### If desired, could fit a GEE that adjusted for baseline to get a better start.

### Inputs:
#### my.data: data to analyze
#### start_theta: starting value for theta
#### start_sigma: starting value for sigma
#### N:
### Output: List of the following:
#### CLWPA_est: estimated effect
#### CLWPA_se: SE of estimate
#### CLWPA_pval: p-value of estimate


CLWPA_fit <- function(my.data, start_theta, start_sigma, N) {

  cml.vertical.adj <- mle2(function(theta,sigma,betabase){loglik_vertical_CML_adj(theta=theta, sigma=sigma,betabase=betabase,
                                                                                  my.data.full=my.data)},start=list(theta=start_theta, sigma=start_sigma,betabase=0))
  #Recreate data frame of all the differences
  #####
  my.data.combos <- data.frame(treated=double(),control=double(),cluster_trt=double(),cluster_ctrl=double())
  for(j in unique(my.data$Period)){
    treated_list.j <- my.data$Summ[my.data$Period == j & my.data$Interv == 1]
    control_list.j <- my.data$Summ[my.data$Period == j & my.data$Interv == 0]
    combos.j <- expand.grid(treated=treated_list.j,control=control_list.j)
    cluster_treated_list.j <- my.data$Cluster[my.data$Period == j & my.data$Interv == 1]
    cluster_control_list.j <- my.data$Cluster[my.data$Period == j & my.data$Interv == 0]
    cluster_combos.j <- expand.grid(cluster_trt=cluster_treated_list.j,cluster_ctrl=cluster_control_list.j)
    my.data.combos <- rbind(my.data.combos,cbind(combos.j,cluster_combos.j))
  }
  my.data.combos$difference <- my.data.combos$treated-my.data.combos$control

  #use cluster ID's in cluster_trt and cluster_ctrl to fill in baseline response values
  #note: merge re-organizes order; that's fine for the CL function
  my.data.combos <- merge(my.data.combos,my.data[my.data$Period == 1,
                                                 c("Cluster","Summ")],by.x="cluster_trt",by.y="Cluster",all.x=TRUE)
  my.data.combos$treated_base <- my.data.combos$Summ
  my.data.combos <- merge(my.data.combos,my.data[my.data$Period == 1,
                                                 c("Cluster","Summ")],by.x="cluster_ctrl",by.y="Cluster",all.x=TRUE)
  my.data.combos$control_base <- my.data.combos$Summ.y

  #difference in baseline
  my.data.combos$difference_base <- my.data.combos$treated_base-my.data.combos$control_base


  #Set up empty matrices to sum
  triangledown_u_sum <- matrix(0,nrow=3,ncol=3)
  u_sq_sum <- matrix(0,nrow=3,ncol=3)



  for(i in 1:nrow(my.data.combos)){
    triangledown_u_sum <- triangledown_u_sum +get_triangledown_u_per_contrast_adj(fitted_params=cml.vertical.adj@coef,
                                                                                  difference=my.data.combos$difference[i],difference_base=my.data.combos$difference_base[i])
  }

  for(i in unique(my.data$Cluster)){
    u.i <- c(0,0,0)#initializing score vector (1x3)
    contrast_involves_i <- my.data.combos$cluster_trt == i | my.data.combos$cluster_ctrl == i
    my.data.combos.i <- my.data.combos[contrast_involves_i,]
    for(k in 1:nrow(my.data.combos.i)){#loop through all the contrasts that involve this cluster
      u.ik <- get_u_per_contrast_adj(fitted_params=cml.vertical.adj@coef,
                                     difference=my.data.combos.i$difference[k],difference_base=my.data.combos.i$difference_base[k])
      u.i <- u.i+u.ik
    }
    u_sq_sum <- u_sq_sum + u.i %*% t(u.i)
  }


  my_H <- - triangledown_u_sum/N
  my_J <- u_sq_sum/N
  my_G <- my_H %*% solve(my_J) %*% my_H
  my_variance_matrix <- solve(my_G)

  return(c(CLWPA_est=cml.vertical.adj@coef["theta"],
           CLWPA_se=((1/N)*my_variance_matrix[1,1])^.5,
           CLWPA_pval=2 * pnorm(abs((cml.vertical.adj@coef["theta"])/(((1/N)*my_variance_matrix[1,1])^.5)),
                                lower.tail = FALSE)))
}
