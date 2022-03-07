#library(tidyverse)
library(dplyr)
library(here)
library(conflicted)
conflicted::conflict_prefer("here", "here")

conflicted::conflict_prefer("filter", "dplyr")

here::here()

job_id <- as.numeric(commandArgs(trailingOnly = T))
task_id <- job_id

iter=1000
b_change = c(0,0.2,0.5,1,2)
error = c(0.2)
b_start = c(-0.5,1,1)
forget_par <- c(0.9958)
obsn_par <- c(1)
updaten_par <- obsn_par

ref_obs_forget <- data.frame(forget_par,obsn_par)

#note, need to expand in reverse order to match random number generation
coef_mat <- expand.grid(run=1:iter,
                        er=error,
                        s2=1:4,
                        b2=b_change,
                        s1=1:4,
                        b1=b_change,
                        obsn_par=obsn_par)

b_change1 <- coef_mat$b1[task_id]
b_change2 <- coef_mat$b2[task_id]
s1 <- coef_mat$s1[task_id]
s2 <- coef_mat$s2[task_id]
er=coef_mat$er[task_id]
it=coef_mat$run[task_id]
forget=ref_obs_forget%>%filter(obsn_par==coef_mat$obsn_par[task_id])%>%pull(forget_par)


#### Run scripts to define the function for each simulation ####
### Linear 1 obs per day
source("DM_sim_linear_1obperday_final_parallelised.R") # DM.sim.linear.1

### Linear N obs per day
source("DM_sim_linear_Nobperday_final.R") # DM.sim.linear.N


#### Run the simulations ####
### Linear 1 obs per day
set.seed(903)


DM.sim.linear.1(it, sampsize=365, valsize = 365,b_change1,b_change2,s1,s2,er,b_start,forget,job_id)





