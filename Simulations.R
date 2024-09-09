#######################################
###### File: Simulations.R ###########
###### Lee Kennedy-Shaffer ############
###### Created 2024/09/09 #############
#######################################

source("A_Const.R")
source("Sigmas.R")
source("Full_Analysis.R")
source("CompEsts.R")

### Parameters for Simulation:
mu <- 0.3
# tau <- 0.06
Alpha1 <- c(-0.009, 0.003, 0.011, -0.019, -0.003, -0.005, -0.012, 
            0.002, 0.005, -0.001, 0.023, 0, 0.018, -0.013)
T1 <- c(0,0.08,0.18,0.29,0.30,0.27,0.20,0.13)
T2 <- c(0,0.02,0.03,0.07,0.13,0.19,0.27,0.30)
ProbT1 <- 1
eta <- 0.01
K <- 100
J <- 8
N <- 14
Theta <- 0

### Framework for simulated data with fixed treatment schedule:
simulate_SWT <- function(NumSims,
                          N, J, 
                          mu, Alpha1, 
                          T1, T2, ProbT1,
                          eta, Theta, K) {
  Sim.Fr <- tibble(Cluster=rep(1:N, each=J),
                    Period=rep(1:J, times=N)) %>%
    dplyr::mutate(Interv=if_else(ceiling(Cluster/2) < Period, 1, 0))
  
  eachsim <- function() {
    Sim.Dat <- Sim.Fr %>%
      dplyr::mutate(FE.g = rep(sample(Alpha1, N, replace=FALSE), each=J),
                    FE.t.Type = rep(sample(c(1,2), size=N, replace=TRUE, prob=c(ProbT1,1-ProbT1)), 
                                    each=J),
                    CPI = rnorm(N*J, mean=0, sd=eta),
                    FE.t=if_else(FE.t.Type==1, T1[Period], T2[Period]),
                    mu.ij=mu+FE.g+FE.t+CPI+Interv*Theta)
    Sim.Dat$Yij <- rbinom(n=N*J, size=rep(K, N*J), prob=pmin(pmax(Sim.Dat$mu.ij,0),1))/rep(K, N*J)
    return(Sim.Dat)
  }
  return(replicate(n=NumSims, as.matrix(eachsim()), simplify="array"))
}

### Generate data:

set.seed(801611)
system.time(a <- simulate_SWT(100, N, J, mu, Alpha1, T1, T2, ProbT1, eta, Theta, K))


