#### simulation function for linear model with 1 observation each day ####

DM.sim.linear.N <- function(iter, sampsize, valsize, filename1, filename2, b_change, error, b_start, obsn, updaten, forget, init_period){
  library(dma)
  library(pROC)
  source("dma_linear_DAJ_source.R")
  library(dplyr)
  #library(plyr)
  # Stop the clock (timer)
  ptm <- proc.time()
  #Create text file to store all the performance results.
  # if(file.exists(paste(filename, ".txt", sep=""))==TRUE){
  #stops from accidentally overwriting a results file.
  #   stop("filename already exists - delete or move file to continue")#.
  # }else{
  filename1 <- here("Results", paste(filename1, ".txt", sep=""))
  OutputNames1 <- c("Scenario1","Scenario2","error","lm_mse", "dm_mse", "dmc_mse", "vc1_mse", "vc2_mse",
                    "b1_lm_mse", "b2_lm_mse", "b3_lm_mse", "b1_dm_mse", "b2_dm_mse", "b3_dm_mse",
                    "b1_dmc_mse", "b2_dmc_mse", "b3_dmc_mse", "b1_vc1_mse", "b2_vc1_mse", "b3_vc1_mse",
                    "b1_vc2_mse"," b2_vc2_mse", "b3_vc2_mse", "b1_true", "b2_true", "b3_true")
  
  filename2 <- here("Results", paste(filename2, ".txt", sep=""))
  
  #Create the output file with the above column headings.
  write(OutputNames1, filename1, append=FALSE, sep="|", ncolumns=length(OutputNames1)) 
  write( c(),filename2, append=FALSE, sep="|", ncolumns=248)
  # }
  
  # setup progress bar
  # pb <- winProgressBar(title = "Iteration Progress bar",
  #                      label = paste("Iterations 0% Completed"),
  #                      min = 0, max = 100, initial = 0)
  
  
  sampsize <- sampsize*obsn
  valsize <- valsize*obsn
  sampsize1 <- sampsize + valsize
  
  # Generate vector detailing which observation number to update the DMs
  update_dev <- seq(round_any(237, updaten, f = ceiling), sampsize, by=updaten)
  update_con <- seq(round_any(237, updaten, f = ceiling), sampsize1, by=updaten)
  if(last(update_dev) != sampsize){update_dev <- c(update_dev,sampsize)}
  if(last(update_con) != sampsize1){update_con <- c(update_con,sampsize1)}
  
  #simulate some data
  #first, static coefficients
  coefmat<-c()
  #coefmat<-cbind(rep(coef[1],sampsize1),rep(coef[2],sampsize1))
  #Add random walk process around static coefficients
  #
  #
  #
  # coef<-c(1.8,0.4)
  # coefmat<-cbind(rep(coef[1],sampsize),rep(coef[2],sampsize))
  
  # Set up change in beta over complete data
  a <- 1 + (valsize/sampsize)
  b <- b_change*a
  lb1 <- length(b)
  
  # generate random number for the intercept and slope of each scenario where b1 is changing over time
  #x1 <- runif(1, 0, 10) # intercept
  #x2 <- runif(1, 0, 10) # slope
  x1 <- b_start[1]
  x2 <- b_start[2]
  x3 <- b_start[3]
  
  # Define dynamic betas for the intercept (1 linear increasing/decreasing, 3 step)
  coefmat1 <- c()
  for(i in 1:lb1){ 
    coefmat1<-cbind(coefmat1,
                    seq(x1,x1+b[i],length.out=sampsize1), # linear increasing
                    c(rep(x1,round(sampsize/2)),rep(x1+b_change[i],round(sampsize1-sampsize/2))) ,
                    c(rep(x1,round(sampsize*.75)),rep(x1+b_change[i],round(sampsize1-sampsize*.75))), 
                    c(rep(x1,round(sampsize*1.083333)),rep(x1+b_change[i],sampsize1-round(sampsize*1.083333)) ) )
    
  }
  coefmat1 <- as.data.frame(coefmat1)
  
  # Define dynamic betas for beta1 (1 linear increasing/decreasing, 3 step)
  coefmat2 <- c()
  for(i in 1:lb1){ 
    coefmat2<-cbind(coefmat2,
                    seq(x2,x2+b[i],length.out=sampsize1), # linear increasing
                    c(rep(x2,round(sampsize/2)),rep(x2+b_change[i],round(sampsize1-sampsize/2))) ,
                    c(rep(x2,round(sampsize*.75)),rep(x2+b_change[i],round(sampsize1-sampsize*.75))), 
                    c(rep(x2,round(sampsize*1.083333)),rep(x2+b_change[i],sampsize1-round(sampsize*1.083333)) ) )
    
  }
  coefmat2 <- as.data.frame(coefmat2)
  coefmat3 <- rep(x3,sampsize1)
  
  
  #lb <- lb1*4 # sets up loop length for all scenarios. For each beta there are 4 scenarios
  
  # coefmat<-cbind(coefmat,seq(1,2.2,length.out=nrow(coefmat)), # linear increaseing
  #                seq(-.75,-1.75,length.out=nrow(coefmat)), # linear decreasing
  #                c(rep(-1.5,nrow(coefmat)/2),rep(-.5,nrow(coefmat)/2)), # step change
  #                c(rep(2.6,nrow(coefmat)/2),rep(1.9,nrow(coefmat)/2))) # step change
  # 
  # coefmat1 <- as.data.frame(coefmat)
  #head(coefmat1)
  
  # Generate time variable
  t1 <- ceiling(sampsize1/updaten)
  t <- sort(rep(seq(1:t1),updaten), decreasing = FALSE)[1:sampsize1]
  t_d <- t[1:sampsize]
  t_v <- t[(sampsize+1):sampsize1]
  #print(length(t))
  
  # mmat1 forces all variables (needed) to be in the model??? setup for dynamic model
  mmat1 <- matrix(c(1,1), nrow=1)
  d <- c()
  npar<-2 # number of predictor variables
  ### loop to generate the coefficient matrix
  for(i in 1:dim(coefmat1)[2]){
    for(k in 1:dim(coefmat2)[2]){
      for(e in seq_along(error)){
        results1 <- vector(mode = "list", length = iter)
        results2 <- vector(mode = "list", length = iter)
        
        for(j in 1:iter){ 
          ## create the data
          d1<-cbind(rnorm(sampsize1,0,1), rnorm(sampsize1,0,1)) # 2 columns of data of length sampsize1
          
          dat_full<-matrix(d1,sampsize1,(npar)) # format data as matrix (needed for DMs)
          dat <- dat_full[1:sampsize,] # subset into training
          dat_val <- dat_full[(sampsize+1):sampsize1,] # subset into validation
          #dat11 <- as.data.frame(dat_full)
          dat1 <- dat_full
          #print(length(dat[,1]))

          # for(i in 1:(lb+1)){
          #   assign(paste0("ydat", i,i), rowSums((cbind(rep(1,nrow(dat_full[,i])),dat_full[,i]))[1:sampsize1,]*coefmat[1:sampsize1,(2*i-1):(2*i)])+rnorm(1, mean = 0, sd = 0.01))
          # }
          intercept_dat <- rep(1,nrow(dat_full))
          # info <- sprintf(paste("Iterations %d%% Completed"), round((j/(iter)*100)))
          # setWinProgressBar(pb, j/((iter))*100, label = info)
          # now create outcome data for each scenario (such that the intercept is the first coefficient used from coefmat)
          #ydat11<-rowSums((cbind(rep(1,nrow(dat_full)),dat_full))[1:sampsize1,]*coefmat[1:sampsize1,1:2])
          assign(paste0("ydat_full"), rowSums(cbind(cbind(intercept_dat,dat_full)[1:sampsize1,]*cbind(coefmat1[,i],coefmat2[,k],coefmat3),
                                                    rnorm(sampsize1, mean = 0, sd = error[e]))))
          
          assign(paste0("ydat"), eval(parse(text=(paste("ydat_full[1:sampsize]",sep="")))))
          
          assign(paste0("ydat_val"), eval(parse(text=(paste("ydat_full[(sampsize+1):sampsize1]",sep="")))))
          
          
          assign(paste0("dma.test"), dma(eval(parse(text=(paste("dat",sep="")))),
                                         eval(parse(text=(paste("ydat",sep="")))),mmat1,lambda=forget,initialperiod=init_period, update_t = update_dev)) # stops at sampsize
          
          
          assign(paste0("dma.test_cont"), dma(eval(parse(text=(paste("dat_full",sep="")))),
                                              eval(parse(text=(paste("ydat_full",sep="")))),mmat1,lambda=forget,initialperiod=init_period, update_t = update_con)) # stops at sampsize
          
          assign(paste0("lm_mod"), lm(eval(parse(text=(paste("ydat ~ dat",sep="")))) ))
          
          assign(paste0("vc_mod1"), lm(eval(parse(text=(paste("ydat ~ dat[,1] + dat[,2] + t_d",sep="")))) ))
          
          assign(paste0("vc_mod2"), lm(eval(parse(text=(paste("ydat ~ dat[,1] + dat[,2] + t_d +
                                                              dat[,1]:t_d + dat[,2]:t_d ",sep="")))) ))
          
          beta <- t(as.data.frame(eval(parse(text=(paste("dma.test$thetahat",sep=""))))))
          
          assign(paste0("dm_b",1), c(rep(beta[sampsize,], valsize)))
          assign(paste0("dm_b",2), c(rep(beta[(2*sampsize),], valsize)))
          assign(paste0("dm_b",3), c(rep(beta[(3*sampsize),], valsize)))
          
          
          beta1 <- t(as.data.frame(eval(parse(text=(paste("dma.test_cont$thetahat",sep=""))))))
          assign(paste0("dm_b_con", 1), c(beta1[sampsize:(sampsize1-1),]) )
          assign(paste0("dm_b_con", 2), c(beta1[(sampsize1+sampsize):(2*sampsize1-1),]))
          assign(paste0("dm_b_con", 3), c(beta1[(2*sampsize1+sampsize):(3*sampsize1-1),]))
          
          assign(paste0("lm_b", 1), eval(parse(text=(paste("lm_mod$coefficients[1]",sep="")))))
          assign(paste0("lm_b", 2), eval(parse(text=(paste("lm_mod$coefficients[2]",sep="")))))
          assign(paste0("lm_b", 3), eval(parse(text=(paste("lm_mod$coefficients[3]",sep="")))))
          
          
          assign(paste0("vc1_b", 1), eval(parse(text=(paste("vc_mod1$coefficients[1] + vc_mod1$coefficients[4]*t_v",sep="")))))
          assign(paste0("vc1_b", 2), eval(parse(text=(paste("vc_mod1$coefficients[2]",sep="")))))
          assign(paste0("vc1_b", 2), c(rep(vc1_b2,valsize)))
          assign(paste0("vc1_b", 3), eval(parse(text=(paste("vc_mod1$coefficients[3]",sep="")))))
          assign(paste0("vc1_b", 3), c(rep(vc1_b3,valsize)))
          
          assign(paste0("vc11_b", 1), eval(parse(text=(paste("vc_mod1$coefficients[1] + vc_mod1$coefficients[4]*rep(sampsize,length(valsize))",sep="")))))
          assign(paste0("vc11_b", 2), eval(parse(text=(paste("vc_mod1$coefficients[2]",sep="")))))
          assign(paste0("vc11_b", 2), c(rep(vc11_b2,valsize)))
          assign(paste0("vc11_b", 3), eval(parse(text=(paste("vc_mod1$coefficients[3]",sep="")))))
          assign(paste0("vc11_b", 3), c(rep(vc11_b3,valsize)))
          
          assign(paste0("vc2_b", 1), eval(parse(text=(paste("vc_mod2$coefficients[1] +", 
                                                            "vc_mod2$coefficients[4]*t_v",sep="")))))
          assign(paste0("vc2_b", 2), eval(parse(text=(paste("vc_mod2$coefficients[2] +", 
                                                            "vc_mod2$coefficients[5]*t_v",sep="")))))
          assign(paste0("vc2_b", 3), eval(parse(text=(paste("vc_mod2$coefficients[3] +", 
                                                            "vc_mod2$coefficients[6]*t_v",sep="")))))
          
          assign(paste0("vc22_b", 1), eval(parse(text=(paste("vc_mod2$coefficients[1] +", 
                                                            "vc_mod2$coefficients[4]*rep(sampsize,length(valsize))",sep="")))))
          assign(paste0("vc22_b", 2), eval(parse(text=(paste("vc_mod2$coefficients[2] +", 
                                                            "vc_mod2$coefficients[5]*rep(sampsize,length(valsize))",sep="")))))
          assign(paste0("vc22_b", 3), eval(parse(text=(paste("vc_mod2$coefficients[3] +", 
                                                            "vc_mod2$coefficients[6]*rep(sampsize,length(valsize))",sep="")))))
          
          # assign(paste0("model_betas"), eval(parse(text=(paste("as.data.frame(cbind(",i,",",k,",error[",e,"],ydat",i,"_",k,"_full[(sampsize+1):sampsize1], dat_full[(sampsize+1):sampsize1,1],dat_full[(sampsize+1):sampsize1,2],",
          #                                                      "coefmat1[(sampsize+1):sampsize1,",i,"],coefmat2[(sampsize+1):sampsize1,",k,"],coefmat3[(sampsize+1):sampsize1],",
          #                                                      "dm_b1 , dm_b2, dm_b3, dm_b_con1,dm_b_con2,dm_b_con3,",
          #                                                      "lm_b1,lm_b2,lm_b3,vc1_b1,vc1_b2,vc1_b3,vc2_b1,vc2_b2,vc2_b3, c((sampsize+1):sampsize1)))",sep="")))))
          # 
          #beta_change <- c( rep(b_change[1], valsize*dim(coefmat1)[2]), rep(b_change[2], valsize*dim(coefmat1)[2]), rep(b_change[3], valsize*dim(coefmat1)[2]),
          #rep(b_change[4], valsize*dim(coefmat1)[2]), rep(b_change[5], valsize*dim(coefmat1)[2]) ) # 
          
          #rep(c(rep("slope",365), rep("step1",365), rep("step2",365), rep("step3",365)),5)
          
          ### Predict outcome in val data for each model
          lm_pred <- lm_b1 + lm_b2*dat_val[,1] + lm_b3*dat_val[,2]
          dm_pred <- dm_b1 + dm_b2*dat_val[,1] + dm_b3*dat_val[,2]
          dmc_pred <- dm_b_con1 + dm_b_con2*dat_val[,1] + dm_b_con3*dat_val[,2]
          vc1_pred <- vc1_b1 + vc1_b2*dat_val[,1] + vc1_b3*dat_val[,2]
          vc2_pred <- vc2_b1 + vc2_b2*dat_val[,1] + vc2_b3*dat_val[,2]
          vc11_pred <- vc11_b1 + vc11_b2*dat_val[,1] + vc11_b3*dat_val[,2]
          vc22_pred <- vc22_b1 + vc22_b2*dat_val[,1] + vc22_b3*dat_val[,2]
          
          ### Calculate MSE (at each timepoint in validation data) for each model
          lm_mse <- (lm_pred - ydat_val)^2
          dm_mse <- (dm_pred - ydat_val)^2
          dmc_mse <- (dmc_pred - ydat_val)^2
          vc1_mse <- (vc1_pred - ydat_val)^2
          vc2_mse <- (vc2_pred - ydat_val)^2
          vc11_mse <- (vc11_pred - ydat_val)^2
          vc22_mse <- (vc22_pred - ydat_val)^2
          
          ### Calculate MSE over full validation data for each model
          m_lm_mse <- mean(lm_mse)
          m_dm_mse <- mean(dm_mse)
          m_dmc_mse <- mean(dmc_mse)
          m_vc1_mse <- mean(vc1_mse)
          m_vc2_mse <- mean(vc2_mse)
          m_vc11_mse <- mean(vc11_mse)
          m_vc22_mse <- mean(vc22_mse)
          
          ### Calculate beta MSE (at each timepoint in validation data) for each model coefmat1[,i],coefmat2[,k]
          b1_lm_mse <- (lm_b1 - coefmat1[(sampsize+1):sampsize1,i])^2
          b2_lm_mse <- (lm_b2 - coefmat2[(sampsize+1):sampsize1,k])^2
          b3_lm_mse <- (lm_b3 - coefmat3[(sampsize+1):sampsize1])^2
          
          b1_dm_mse <- (dm_b1 - coefmat1[(sampsize+1):sampsize1,i])^2
          b2_dm_mse <- (dm_b2 - coefmat2[(sampsize+1):sampsize1,k])^2
          b3_dm_mse <- (dm_b3 - coefmat3[(sampsize+1):sampsize1])^2
          
          b1_dmc_mse <- (dm_b_con1 - coefmat1[(sampsize+1):sampsize1,i])^2
          b2_dmc_mse <- (dm_b_con2 - coefmat2[(sampsize+1):sampsize1,k])^2
          b3_dmc_mse <- (dm_b_con3 - coefmat3[(sampsize+1):sampsize1])^2
          
          b1_vc1_mse <- (vc1_b1 - coefmat1[(sampsize+1):sampsize1,i])^2
          b2_vc1_mse <- (vc1_b2 - coefmat2[(sampsize+1):sampsize1,i])^2
          b3_vc1_mse <- (vc1_b3 - coefmat3[(sampsize+1):sampsize1])^2
          
          b1_vc2_mse <- (vc2_b1 - coefmat1[(sampsize+1):sampsize1,i])^2
          b2_vc2_mse <- (vc2_b2 - coefmat2[(sampsize+1):sampsize1,k])^2
          b3_vc2_mse <- (vc2_b3 - coefmat3[(sampsize+1):sampsize1])^2
          
          
          ### Calculate beta MSE over full validation data for each model
          m_b1_lm_mse <- mean(b1_lm_mse)
          m_b2_lm_mse <- mean(b2_lm_mse)
          m_b3_lm_mse <- mean(b3_lm_mse)
          
          m_b1_dm_mse <- mean(b1_dm_mse)
          m_b2_dm_mse <- mean(b2_dm_mse)
          m_b3_dm_mse <- mean(b3_dm_mse)
          
          m_b1_dmc_mse <- mean(b1_dmc_mse)
          m_b2_dmc_mse <- mean(b2_dmc_mse)
          m_b3_dmc_mse <- mean(b3_dmc_mse)
          
          m_b1_vc1_mse <- mean(b1_vc1_mse)
          m_b2_vc1_mse <- mean(b2_vc1_mse)
          m_b3_vc1_mse <- mean(b3_vc1_mse)
          
          m_b1_vc2_mse <- mean(b1_vc2_mse)
          m_b2_vc2_mse <- mean(b2_vc2_mse)
          m_b3_vc2_mse <- mean(b3_vc2_mse)
          
          ### Calibration slope
          cal_s_lm_mse <- lm(ydat_val ~ lm_pred) 
          cal_s_dm_mse <- lm(ydat_val ~ dm_pred)
          cal_s_dmc_mse <- lm(ydat_val ~ dmc_pred)
          cal_s_vc1_mse <- lm(ydat_val ~ vc1_pred)
          cal_s_vc2_mse <- lm(ydat_val ~ vc2_pred)
          cal_s_vc11_mse <- lm(ydat_val ~ vc11_pred)
          cal_s_vc22_mse <- lm(ydat_val ~ vc22_pred)
          
          ### Calibration intercept
          cal_i_lm_mse <- lm(ydat_val ~ 1, offset = lm_pred) 
          cal_i_dm_mse <- lm(ydat_val ~ 1, offset = dm_pred)
          cal_i_dmc_mse <- lm(ydat_val ~ 1, offset = dmc_pred)
          cal_i_vc1_mse <- lm(ydat_val ~ 1, offset = vc1_pred)
          cal_i_vc2_mse <- lm(ydat_val ~ 1, offset = vc2_pred)
          cal_i_vc11_mse <- lm(ydat_val ~ 1, offset = vc11_pred)
          cal_i_vc22_mse <- lm(ydat_val ~ 1, offset = vc22_pred)
          
          ### disc
          disc_lm_mse <- roc(ydat_val ~ lm_pred, quiet = TRUE) 
          disc_dm_mse <- roc(ydat_val ~ dm_pred, quiet = TRUE)
          disc_dmc_mse <- roc(ydat_val ~ dmc_pred, quiet = TRUE)
          disc_vc1_mse <- roc(ydat_val ~ vc1_pred, quiet = TRUE)
          disc_vc2_mse <- roc(ydat_val ~ vc2_pred, quiet = TRUE)
          disc_vc11_mse <- roc(ydat_val ~ vc11_pred, quiet = TRUE)
          disc_vc22_mse <- roc(ydat_val ~ vc22_pred, quiet = TRUE)
          
          
          ###
          results1[[j]] <- cbind(lm_mse, dm_mse, dmc_mse, vc1_mse, vc2_mse,
                                 b1_lm_mse, b2_lm_mse, b3_lm_mse,
                                 b1_dm_mse, b2_dm_mse, b3_dm_mse,
                                 b1_dmc_mse, b2_dmc_mse, b3_dmc_mse,
                                 b1_vc1_mse, b2_vc1_mse, b3_vc1_mse,
                                 b1_vc2_mse, b2_vc2_mse, b3_vc2_mse)
          row.names(results1[[j]])<- NULL ; colnames(results1[[j]])<- NULL
          #results1[[j]] <- as.matrix(results1[[j]])
          
          results2[[j]] <- c(as.numeric(cal_s_lm_mse$coefficients[2]), as.numeric(cal_s_dm_mse$coefficients[2]), as.numeric(cal_s_dmc_mse$coefficients[2]), as.numeric(cal_s_vc1_mse$coefficients[2]), as.numeric(cal_s_vc2_mse$coefficients[2]), as.numeric(cal_s_vc11_mse$coefficients[2]), as.numeric(cal_s_vc22_mse$coefficients[2]),
                             as.numeric(cal_i_lm_mse$coefficients[1]), as.numeric(cal_i_dm_mse$coefficients[1]), as.numeric(cal_i_dmc_mse$coefficients[1]), as.numeric(cal_i_vc1_mse$coefficients[1]), as.numeric(cal_i_vc2_mse$coefficients[1]), as.numeric(cal_i_vc11_mse$coefficients[1]), as.numeric(cal_i_vc22_mse$coefficients[1]),
                             disc_lm_mse$auc, disc_dm_mse$auc, disc_dmc_mse$auc, disc_vc1_mse$auc, disc_vc2_mse$auc, disc_vc11_mse$auc, disc_vc22_mse$auc,
                             m_lm_mse, m_dm_mse, m_dmc_mse, m_vc1_mse, m_vc2_mse, m_vc11_mse, m_vc22_mse)
          
          #results[j] <- cbind(j, model_betas)
          #write.table(all_betas, filename, append=TRUE, sep="|", row.names=FALSE, col.names=FALSE) 
        }
        ### Take results1 and average over all iterations
        #r1 = laply(results1, as.matrix)
        r1 = aaply(laply(results1, as.matrix), c(2, 3), mean)
        r11 <- cbind(i, k, error[e], r1, coefmat1[(sampsize+1):sampsize1,i],
                     coefmat2[(sampsize+1):sampsize1,k], coefmat3[(sampsize+1):sampsize1])
        
        ### Take results2 and average over iterations
        # average over column and then create row of data including scenario and output to results
        r2 <- do.call(rbind,results2)
        # Obtain summary info across the 1000 iteratons - this generates vector of length 213 (3 + (30*7))
        r22 <- matrix(c(i, k, error[e],
                        apply(r2, 2, FUN = mean),
                        apply(r2, 2, FUN = sd),
                        apply(r2 ,2, FUN = median),
                        apply(r2 ,2, FUN = quantile, probs=c(0.025,0.1,0.9,0.975)) ), nrow=1)
        
        # Output as row to file ...
        #write.table(r11, filename1, append=TRUE, sep="|", row.names=FALSE, col.names=FALSE) 
        write.table(r22, filename2, append=TRUE, sep="|", row.names=FALSE, col.names=FALSE) 
        
        # 
      }
      print(c(i,k))
      
    }
    
    #print(all_betas[1,])
    #print(dim(model_betas))
    #print(length(lm_b1))
    
    #print(length(vc1_b1))
    #print(length(vc2_b2))
    # all_betas <- cbind(j, model_betas)
    
  }
  #close(pb)
  #print(dim(r1))
  #print(r1)
  #print(dim(r11))
  #print(r22)
  
  # Stop the clock
  proc.time() - ptm
}

