#### simulation function for linear model with 1 observation each day ####

DM.sim.linear.1 <- function(iter, sampsize=365, valsize = 365,b_change1,b_change2,s1,s2,er,b_start,forget,job_id){
  library(dma)
  library(pROC)
  #library(plyr)
  # Stop the clock (timer)
  ptm <- proc.time()

  # setup progress bar
  # pb <- winProgressBar(title = "Iteration Progress bar",
  #                      label = paste("Iterations 0% Completed"),
  #                      min = 0, max = 100, initial = 0)
  
  #make sure random number generation starts from correct point
  if(task_id>0) {rnorm(3*task_id*sampsize1)}

    
    sampsize1 <- sampsize + valsize
    #simulate some data
    #first, static coefficients
    coefmat<-c()
      
    # Set up change in beta over complete data
    a <- 1 + (valsize/sampsize)
    b1 <- b_change1*a
    b2 <- b_change2*a
    
    x1 <- b_start[1]
    x2 <- b_start[2]
    x3 <- b_start[3]
    
    # Define dynamic betas for the intercept (1 linear increasing/decreasing, 3 step)
      if(s1==1){coefmat1 <- seq(x1,x1+b1,length.out=sampsize1)}
      if(s1==2){coefmat1 <- c(rep(x1,round(sampsize/2)),rep(x1+b_change1,round(sampsize1-sampsize/2)))}
      if(s1==3){coefmat1 <- c(rep(x1,round(sampsize*.75)),rep(x1+b_change1,round(sampsize1-sampsize*.75)))} 
      if(s1==4){coefmat1 <- c(rep(x1,round(sampsize*1.083333)),rep(x1+b_change1,sampsize1-round(sampsize*1.083333)))}
    
    coefmat1 <- as.data.frame(coefmat1)
    
    # Define dynamic betas for beta1 (1 linear increasing/decreasing, 3 step)
    if(s2==1){coefmat2 <- seq(x2,x2+b2,length.out=sampsize1)}
    if(s2==2){coefmat2 <- c(rep(x2,round(sampsize/2)),rep(x2+b_change2,round(sampsize1-sampsize/2)))}
    if(s2==3){coefmat2 <- c(rep(x2,round(sampsize*.75)),rep(x2+b_change2,round(sampsize1-sampsize*.75)))} 
    if(s2==4){coefmat2 <- c(rep(x2,round(sampsize*1.083333)),rep(x2+b_change2,sampsize1-round(sampsize*1.083333)))}
    
    coefmat2 <- as.data.frame(coefmat2)
    
    coefmat3 <- rep(x3,sampsize1)

    
    # mmat1 forces all variables (needed) to be in the model??? setup for dynamic model
    mmat1 <- matrix(c(1,1), nrow=1)
    d <- c()
    npar<-2 # number of predictor variables
   
            ## create the data
            d1<-cbind(rnorm(sampsize1,0,1), rnorm(sampsize1,0,1)) # 2 columns of data of length sampsize1
            
            dat_full<-matrix(d1,sampsize1,(npar)) # format data as matrix (needed for DMs)
            dat <- dat_full[1:sampsize,] # subset into training
            dat_val <- dat_full[(sampsize+1):sampsize1,] # subset into validation
            #dat11 <- as.data.frame(dat_full)
            dat1 <- dat_full
            
            ## Generate time variable
            t <- 1:sampsize
            
            # for(i in 1:(lb+1)){
            #   assign(paste0("ydat", i,i), rowSums((cbind(rep(1,nrow(dat_full[,i])),dat_full[,i]))[1:sampsize1,]*coefmat[1:sampsize1,(2*i-1):(2*i)])+rnorm(1, mean = 0, sd = 0.01))
            # }
            intercept_dat <- rep(1,nrow(dat_full))
            # info <- sprintf(paste("Iterations %d%% Completed"), round((j/(iter)*100)))
            # setWinProgressBar(pb, j/((iter))*100, label = info)
          # now create outcome data for each scenario (such that the intercept is the first coefficient used from coefmat)
          #ydat11<-rowSums((cbind(rep(1,nrow(dat_full)),dat_full))[1:sampsize1,]*coefmat[1:sampsize1,1:2])
            assign(paste0("ydat_full"), rowSums(cbind(cbind(intercept_dat,dat_full)[1:sampsize1,]*cbind(coefmat1,coefmat2,coefmat3),
                                                    rnorm(sampsize1, mean = 0, sd = er))))
          
          assign(paste0("ydat"), eval(parse(text=(paste("ydat_full[1:sampsize]",sep="")))))
          
          assign(paste0("ydat_val"), eval(parse(text=(paste("ydat_full[(sampsize+1):sampsize1]",sep="")))))
      
          
          assign(paste0("dma.test"), dma(eval(parse(text=(paste("dat",sep="")))),
                                            eval(parse(text=(paste("ydat",sep="")))),mmat1,lambda=forget,initialperiod=1)) # stops at sampsize
          

          assign(paste0("dma.test_cont"), dma(eval(parse(text=(paste("dat_full",sep="")))),
                                                       eval(parse(text=(paste("ydat_full",sep="")))),mmat1,lambda=forget,initialperiod=1)) # stops at sampsize
          
          assign(paste0("lm_mod"), lm(eval(parse(text=(paste("ydat ~ dat",sep="")))) ))
          
          assign(paste0("vc_mod1"), lm(eval(parse(text=(paste("ydat ~ dat[,1] + dat[,2] + t",sep="")))) ))
          
          assign(paste0("vc_mod2"), lm(eval(parse(text=(paste("ydat ~ dat[,1] + dat[,2] + t +
                                                              dat[,1]:t + dat[,2]:t ",sep="")))) ))

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
          
          
          assign(paste0("vc1_b", 1), eval(parse(text=(paste("vc_mod1$coefficients[1] + vc_mod1$coefficients[4]*c((sampsize+1):sampsize1)",sep="")))))
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
                                                            "vc_mod2$coefficients[4]*c((sampsize+1):sampsize1)",sep="")))))
          assign(paste0("vc2_b", 2), eval(parse(text=(paste("vc_mod2$coefficients[2] +", 
                                                            "vc_mod2$coefficients[5]*c((sampsize+1):sampsize1)",sep="")))))
          assign(paste0("vc2_b", 3), eval(parse(text=(paste("vc_mod2$coefficients[3] +", 
                                                            "vc_mod2$coefficients[6]*c((sampsize+1):sampsize1)",sep="")))))
          
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
          b1_lm_mse <- (lm_b1 - coefmat1[(sampsize+1):sampsize1,1])^2
          b2_lm_mse <- (lm_b2 - coefmat2[(sampsize+1):sampsize1,1])^2
          b3_lm_mse <- (lm_b3 - coefmat3[(sampsize+1):sampsize1])^2
          
          b1_dm_mse <- (dm_b1 - coefmat1[(sampsize+1):sampsize1,1])^2
          b2_dm_mse <- (dm_b2 - coefmat2[(sampsize+1):sampsize1,1])^2
          b3_dm_mse <- (dm_b3 - coefmat3[(sampsize+1):sampsize1])^2
          
          b1_dmc_mse <- (dm_b_con1 - coefmat1[(sampsize+1):sampsize1,1])^2
          b2_dmc_mse <- (dm_b_con2 - coefmat2[(sampsize+1):sampsize1,1])^2
          b3_dmc_mse <- (dm_b_con3 - coefmat3[(sampsize+1):sampsize1])^2
          
          b1_vc1_mse <- (vc1_b1 - coefmat1[(sampsize+1):sampsize1,1])^2
          b2_vc1_mse <- (vc1_b2 - coefmat2[(sampsize+1):sampsize1,1])^2
          b3_vc1_mse <- (vc1_b3 - coefmat3[(sampsize+1):sampsize1])^2
          
          b1_vc2_mse <- (vc2_b1 - coefmat1[(sampsize+1):sampsize1,1])^2
          b2_vc2_mse <- (vc2_b2 - coefmat2[(sampsize+1):sampsize1,1])^2
          b3_vc2_mse <- (vc2_b3 - coefmat3[(sampsize+1):sampsize1])^2
          
          b1_vc11_mse <- (vc11_b1 - coefmat1[(sampsize+1):sampsize1,1])^2
          b2_vc11_mse <- (vc11_b2 - coefmat2[(sampsize+1):sampsize1,1])^2
          b3_vc11_mse <- (vc11_b3 - coefmat3[(sampsize+1):sampsize1])^2
          
          b1_vc22_mse <- (vc22_b1 - coefmat1[(sampsize+1):sampsize1,1])^2
          b2_vc22_mse <- (vc22_b2 - coefmat2[(sampsize+1):sampsize1,1])^2
          b3_vc22_mse <- (vc22_b3 - coefmat3[(sampsize+1):sampsize1])^2
          
          ### Calculate beta MSE over full validation data for each model
          # m_b1_lm_mse <- mean(b1_lm_mse)
          # m_b2_lm_mse <- mean(b2_lm_mse)
          # m_b3_lm_mse <- mean(b3_lm_mse)
          # 
          # m_b1_dm_mse <- mean(b1_dm_mse)
          # m_b2_dm_mse <- mean(b2_dm_mse)
          # m_b3_dm_mse <- mean(b3_dm_mse)
          # 
          # m_b1_dmc_mse <- mean(b1_dmc_mse)
          # m_b2_dmc_mse <- mean(b2_dmc_mse)
          # m_b3_dmc_mse <- mean(b3_dmc_mse)
          # 
          # m_b1_vc1_mse <- mean(b1_vc1_mse)
          # m_b2_vc1_mse <- mean(b2_vc1_mse)
          # m_b3_vc1_mse <- mean(b3_vc1_mse)
          # 
          # m_b1_vc2_mse <- mean(b1_vc2_mse)
          # m_b2_vc2_mse <- mean(b2_vc2_mse)
          # m_b3_vc2_mse <- mean(b3_vc2_mse)
          # 
          # m_b1_vc11_mse <- mean(b1_vc11_mse)
          # m_b2_vc11_mse <- mean(b2_vc11_mse)
          # m_b3_vc11_mse <- mean(b3_vc11_mse)
          # 
          # m_b1_vc22_mse <- mean(b1_vc22_mse)
          # m_b2_vc22_mse <- mean(b2_vc22_mse)
          # m_b3_vc22_mse <- mean(b3_vc22_mse)
          
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
          
          # disc_lm_mse <- disc_lm_mse1$auc
          # disc_dm_mse <- disc_dm_mse1$auc
          # disc_dmc_mse <- disc_dmc_mse1$auc
          # disc_vc1_mse <- disc_vc1_mse1$auc
          # disc_vc2_mse <- disc_vc2_mse1$auc
          
          
          ###
          results1 <- cbind(lm_mse, dm_mse, dmc_mse, vc1_mse, vc2_mse, vc11_mse, vc22_mse,
                                 b1_lm_mse, b2_lm_mse, b3_lm_mse,
                                 b1_dm_mse, b2_dm_mse, b3_dm_mse,
                                 b1_dmc_mse, b2_dmc_mse, b3_dmc_mse,
                                 b1_vc1_mse, b2_vc1_mse, b3_vc1_mse,
                                 b1_vc2_mse, b2_vc2_mse, b3_vc2_mse,
                            b1_vc11_mse, b2_vc11_mse, b3_vc11_mse,
                            b1_vc22_mse, b2_vc22_mse, b3_vc22_mse)
          row.names(results1)<- NULL ; colnames(results1)<- NULL
          #results1[[j]] <- as.matrix(results1[[j]])
          
          results2 <- c(as.numeric(cal_s_lm_mse$coefficients[2]), as.numeric(cal_s_dm_mse$coefficients[2]), as.numeric(cal_s_dmc_mse$coefficients[2]), as.numeric(cal_s_vc1_mse$coefficients[2]), as.numeric(cal_s_vc2_mse$coefficients[2]), as.numeric(cal_s_vc11_mse$coefficients[2]), as.numeric(cal_s_vc22_mse$coefficients[2]),
                             as.numeric(cal_i_lm_mse$coefficients[1]), as.numeric(cal_i_dm_mse$coefficients[1]), as.numeric(cal_i_dmc_mse$coefficients[1]), as.numeric(cal_i_vc1_mse$coefficients[1]), as.numeric(cal_i_vc2_mse$coefficients[1]), as.numeric(cal_i_vc11_mse$coefficients[1]), as.numeric(cal_i_vc22_mse$coefficients[1]),
                             disc_lm_mse$auc, disc_dm_mse$auc, disc_dmc_mse$auc, disc_vc1_mse$auc, disc_vc2_mse$auc, disc_vc11_mse$auc, disc_vc22_mse$auc,
                             m_lm_mse, m_dm_mse, m_dmc_mse, m_vc1_mse, m_vc2_mse, m_vc11_mse, m_vc22_mse)
          
          saveRDS(results1,paste0("results1_b1_",b_change1,"_scen_",s1,"_b2_",b_change2,"_scen_",s2,"_sd_",e,"_iter_",iter,".RDS"))
          saveRDS(results2,paste0("results2_b1_",b_change1,"_scen_",s1,"_b2_",b_change2,"_scen_",s2,"_sd_",e,"_iter_",iter,".RDS"))
          
          #results[j] <- cbind(j, model_betas)
          #write.table(all_betas, filename, append=TRUE, sep="|", row.names=FALSE, col.names=FALSE) 
          }
          ### Take results1 and average over all iterations
          #r1 = laply(results1, as.matrix)
          # r1 = aaply(laply(results1, as.matrix), c(2, 3), mean)
          # r11 <- cbind(i, k, error[e], r1, coefmat1[(sampsize+1):sampsize1,i],
          #          coefmat2[(sampsize+1):sampsize1,k], coefmat3[(sampsize+1):sampsize1])
          # 
          # ### Take results2 and average over iterations
          # # average over column and then create row of data including scenario and output to results
          # r2 <- do.call(rbind,results2)
          # # Obtain summary info across the 1000 iteratons - this generates vector of length 248 (3 + (35*7))
          # r22 <- matrix(c(i, k, error[e],
          # apply(r2, 2, FUN = mean),
          # apply(r2, 2, FUN = sd),
          # apply(r2 ,2, FUN = median),
          # apply(r2 ,2, FUN = quantile, probs=c(0.025,0.1,0.9,0.975)) ), nrow=1)
          # 
          # # Output as row to file ...
          # #write.table(r11, filename1, append=TRUE, sep="|", row.names=FALSE, col.names=FALSE) 
          # write.table(r22, filename2, append=TRUE, sep="|", row.names=FALSE, col.names=FALSE) 
          
          # 
     
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

