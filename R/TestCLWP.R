## Test Sims with CLWP ##
source("R/Simulations.R")

index <- 1
Sim.Dat <- Sim_Data(Sim.Fr=Sim_Frame(Param_Set$N[index],Param_Set$J[index]),
                    mu=Param_Set$mu[index],
                    Alpha1=Alpha1,
                    T1=T1,
                    T2=T2,
                    ProbT1=Param_Set$ProbT1[index],
                    sig_nu=Param_Set$sig_nu[index],
                    sig_e=Param_Set$sig_e[index],
                    m=Param_Set$m[index],
                    ThetaType=Theta_Set[[index]]$Type,
                    ThetaDF=Theta_Set[[index]]$ThetaDF)
