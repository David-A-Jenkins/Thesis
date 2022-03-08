#### SPC simulation ####

### First load the BCIS data 
library(here)
library(dplyr)
df <- readRDS(here("BCIS_updated_cleaned_subset.Rda"))

### Generate LP using coefficients from Kate's paper
lp <- -6.089 +0.071*df$age + 0.114*df$sex + 0.524*df$diabetes +
  0.158*df$mi + 0.997*df$renal1 + 1.128*df$renal2 + 0.43*df$ce +
  1.004*df$indication1 + 2.114*df$indication2 + 2.295*df$indication3 + 2.531*df$indication4 +
  3.817*df$cs - 0.026*df$age_shock - 0.951*df$indication_shock1 - 
  1.226*df$indication_shock2 - 1.203*df$indication_shock3 - 1.438*df$indication_shock4 - 0.016*df$age_diabetes


trig <- c()

drift <- seq(0,2,0.1) # set up miscalibration-in-the-large

for(b in seq_along(drift)){ # loop over intercept drift (miscalibration-in-thelarge) scenarios
  print(drift[b])
  
  for(i in 1:1000){ # loop over iterations
    set.seed(i)
    lp_mod <- sample(lp, replace=TRUE) # Samples from the data with replacement
    
    assign(paste0("lp_true"), lp_mod+drift[b]) # Generates the 'true' LP by adding the drift to the intercept
    
    pr_true <- exp(lp_true)/(1+exp(lp_true)) # generate 'true' probability 
    pr_mod <- exp(lp_mod)/(1+exp(lp_mod)) # generate model probability 
    
    var_pr <- pr_mod*(1-pr_mod) # Variance
    var_pr_sum <- cumsum(var_pr) # Variance is addative and not SD
    sd_pr_sum <- sqrt(var_pr_sum) # sqrt the sum of variance to get SD
    lci3 <- -3*sd_pr_sum # Generate lower limit
    uci3 <- 3*sd_pr_sum # generate upper limit
    lci4 <- -4*sd_pr_sum # Generate lower limit
    uci4 <- 4*sd_pr_sum # generate upper limit
    
    ## Generate true outcome
    dat_true_p <- rbinom(length(pr_true),1, pr_true)
    
    ## Generate O-E
    OE <- dat_true_p - pr_mod
    
    ## Generate control limits and cumsum of O-E
    spc_OE <- cumsum(OE)
    
    spc_OE <- as.data.frame(cbind(lci3, uci3, lci4, uci4, spc_OE)) # combine data
    
    
    a <- min(which(spc_OE$spc_OE<lci3 & cumsum(pr_mod)>5))
    a1 <- min(which(spc_OE$spc_OE>uci3 & cumsum(pr_mod)>5))
    a2 <- min(which(spc_OE$spc_OE<lci4 & cumsum(pr_mod)>5))
    a3 <- min(which(spc_OE$spc_OE>uci4 & cumsum(pr_mod)>5))
    
    t1 <- c(drift[b], i, min(a,a1), min(a2,a3))
    
    trig <- rbind(trig, t1)
    
  }
}
trig <- as.data.frame(trig)
trig <- trig %>% rename(drift=V1, iter=V2, trig3=V3, trig4=V4) # trig3 is the observation number of alert (if there was an alert)

saveRDS(trig, file = "trig_int_noreset_exp5.rds") # saves file

#trig <- readRDS(file = "trig_null_test.rds")


#
#
#


#


#
#
#
#
