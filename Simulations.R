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
tau <- 0.06
T1 <- c(0,0.08,0.18,0.29,0.30,0.27,0.20,0.13)
T2 <- c(0,0.02,0.03,0.07,0.13,0.19,0.27,0.30)
eta <- 0.01
K <- 100
J <- 8
N <- 14