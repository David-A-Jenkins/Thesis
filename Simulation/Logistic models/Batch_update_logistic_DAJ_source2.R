# ----------------------------------------------------

## Updated: logistic.dma.default

logistic.dma.batch <-
  function (x, y, models.which, lambda = 0.99, alpha = 0.99, autotune = TRUE, 
            initmodelprobs = NULL, initialsamp = NULL, bite.size) 
  {
    require(mnormt)
    require(MASS)
    K <- nrow(models.which) #number of models considered
    baytah <- array(dim = c(nrow(models.which), nrow(x), (ncol(x) + 
                                                            1)))
    varbaytah <- array(dim = c(nrow(models.which), nrow(x), (ncol(x) + 
                                                               1)))
    laplacemodel <- array(dim = c(nrow(models.which), nrow(x)))
    yhatmodel <- array(dim = c(nrow(models.which), nrow(x)))
    if (is.null(initialsamp)) {
      initialsamp <- round(nrow(x)/10, 0) #if no amount of initial samples are given use first 10%
    }
    for (mm in 1:nrow(models.which)) {
     # pb <- winProgressBar(title="Progress bar", label=paste("0% done on model", mm), min=0, max=100, initial=0)
      xdat <- x[, c(which(models.which[mm, ] == 1))] #all covariates which are included in model mm
      
      xdat <- cbind(1, xdat) #adds an intercept to the design matrix
      d <- dim(xdat)[2] #number of covaraites in the model
      if (is.matrix(autotune)) {
        tune.mat <- autotune #tests if the autotune variable is a matrix or not
      }
      if (!is.matrix(autotune)) { #autotune decides if we chooce between forgetting and no forgetting at each time point. 
        if (autotune == TRUE) {
          tune.mat <- tunemat.fn(lambda, 1, d)
        }
        if (autotune == FALSE) {            
          tune.mat <- matrix(rep(lambda, d), nrow = 1, ncol = d)
        }
      }
      init.temp <- dlogr.init(xdat[1:initialsamp, ], y[1:initialsamp]) #dlogr.init function fits a log reg model with passed data 
      betahat.tm1 <- init.temp$BetaHat #this is the estimated coefficients from the glm with the initial sample of data
      baytah[mm, 1:initialsamp, c(1, 1 + which(models.which[mm,] == 1))] <- matrix(rep(init.temp$BetaHat, initialsamp), 
                                                                                initialsamp, length(init.temp$BetaHat), 
                                                                                byrow = TRUE) #this puts the coefficient estimates
      #from the initial sample into the appropriate places in the array baytah
      varbaytah[mm, 1:initialsamp, c(1, 1 + which(models.which[mm, 
                                                               ] == 1))] <- matrix(rep(diag(init.temp$VarBetaHat), 
                                                                                       initialsamp), initialsamp, length(init.temp$BetaHat), 
                                                                                   byrow = TRUE) #this puts the var estimates
      #from the initial sample into the appropriate places in the array varbaytah
      varbetahat.tm1 <- init.temp$VarBetaHat
      laplacemodel[mm, 1:initialsamp] <- rep(0, initialsamp)
      yhatmodel[mm, 1:initialsamp] <- exp(xdat[1:initialsamp,] %*% init.temp$BetaHat)/
                                                                    (1 + exp(xdat[1:initialsamp,] %*% init.temp$BetaHat)) #puts the
      #estimated risks into the first inital sample place for model mm
      
      for (i in 1:(length(bite.size)-1)) { #bite size is a vector outlining how many observations are in each update batch; the first
        #position gives how many are used in the initial sample- i.e. bite size gives the posiitons in the data for which all the 
        #different updates take place for instance c(12,24,30,50,77) would use data in rows 1:12 for initial then use 13:24 for the 
        #first update; 25:30 for the second update etc.
        #print(i)
        x.t <- (xdat[(bite.size[i] + 1):(bite.size[i+1]), ])
        y.t <- y[(bite.size[i] + 1):(bite.size[i+1])]
        step.tmp <- dlogr.step(x.t, y.t, betahat.tm1, varbetahat.tm1, 
                               tune.mat)
        
        baytah[mm, (bite.size[i] + 1):(bite.size[i+1]), c(1, 1 + which(models.which[mm,] == 1))] <- 
                                                                      rep(step.tmp$betahat.t, each = (bite.size[i+1]-bite.size[i]))
        varbaytah[mm, (bite.size[i] + 1):(bite.size[i+1]), c(1, 1 + which(models.which[mm,] == 1))] <- 
                                                              rep(diag(step.tmp$varbetahat.t), each = (bite.size[i+1]-bite.size[i]))
        
        laplacemodel[mm, (bite.size[i] + 1):(bite.size[i+1])] <- step.tmp$laplace.t
        
        yhatmodel[mm, (bite.size[i] + 1):(bite.size[i+1])] <- exp(x.t %*% step.tmp$betahat.t)/(1 + 
                                                                                                 exp(x.t %*% step.tmp$betahat.t))
        betahat.tm1 <- step.tmp$betahat.t # replaces the current estimate of beta as the previous values for use next for loop
        varbetahat.tm1 <- step.tmp$varbetahat.t
        lambda = rbind(lambda, step.tmp$lambda)
        
       # info <- sprintf(paste("%d%% done on model", mm), round((i/(length(bite.size)-1))*100))
       # setWinProgressBar(pb, i/((length(bite.size)-1))*100, label=info)
      }
     # close(pb)
     # print(paste("Finished processing model", mm))
    }
    if (sum(is.na(laplacemodel)) > 0 | sum(laplacemodel == Inf) > 
          0 | sum(laplacemodel == -Inf) > 0) {
      print("Warning: At least one laplace approximation is not well behaved. This will likely lead to issues with posterior model probabilities.\n This is likely a computation issue")
    }
    pi.t.t <- array(dim = c(K, nrow(x)))#, dimnames = list(,c("model", "time")))
    pi.t.tm1 <- array(dim = c(K, nrow(x)))#, dimnames = list(c("model", "time")))
    omega.tk <- array(dim = c(K, nrow(x)))#, dimnames = list(c("model", "time")))
    if (is.null(initmodelprobs)) {
      pi.t.t[, 1] <- 1/K
      pi.t.tm1[, 1] <- 1/K
      omega.tk[, 1] <- 1/K
    }
    else {
      pi.t.t[, 1] <- initmodelprobs
      pi.t.tm1[, 1] <- initmodelprobs
      omega.tk[, 1] <- initmodelprobs
    }
    alpha.vec <- rep(NA, nrow(x))
    alpha.vec[1] <- alpha
    alpha.noforget = 1
    for (t in 2:length(alpha.vec)) {
      if (t > 2) {
        rm(alpha.tmp)
      }
      pred.forget <- sum(exp(laplacemodel[, t]) * ((pi.t.t[, 
                                                           t - 1]^alpha)/sum(pi.t.t[, t - 1]^alpha)))
      pred.noforget <- sum(exp(laplacemodel[, t]) * ((pi.t.t[, 
                                                             t - 1]^alpha.noforget)/sum(pi.t.t[, t - 1]^alpha.noforget)))
      alpha.vec[t] <- ifelse(pred.forget >= pred.noforget, 
                             alpha, alpha.noforget)
      alpha.tmp <- alpha.vec[t]
      for (k in 1:K) {
        pi.t.tm1[k, t] <- (pi.t.t[k, t - 1]^alpha.tmp)/sum(pi.t.t[, 
                                                                  t - 1]^alpha.tmp)
        omega.tk[k, t] <- (pi.t.tm1[k, t] * exp(laplacemodel[k, 
                                                             t]))
      }
      for (k in 1:K) {
        pi.t.t[k, t] <- omega.tk[k, t]/sum(omega.tk[, t])
        pi.t.t[k, t] <- ifelse(pi.t.t[k, t] > 0, pi.t.t[k, 
                                                        t], 0.001)
        pi.t.t[k, t] <- ifelse(pi.t.t[k, t] < 1, pi.t.t[k, 
                                                        t], 0.999)
      }
    }
    yhatdma = colSums(yhatmodel * pi.t.t)
    est <- (list(x = x, y = y, models = models.which, lambda = lambda, 
                 alpha = alpha, pmp = pi.t.t, alpha.used = alpha.vec, 
                 theta = baytah, vartheta = varbaytah, yhatdma = yhatdma, 
                 yhatmodel = yhatmodel))
    est$fitted.values <- yhatdma
    est$residuals <- y - yhatdma
    class(est) <- "logistic.dma"
    est
  }


#------------------------------------------------------------
## Updated: tunemat.fn  

tunemat.fn <- 
  function (l.forget, l.noforget, d) #sclaer of forgetting, no forgetting and number of covariates respectively
  {
    bin <- 1:(2^d) #2^d is all possible combinations of analysing 'no forgetting' and 'some forgetting' at each time point.
    tune.mat.temp <- sapply(1:d, function(i) { #puts the values 1:d into the written function
      r <- as.logical(bin %% 2)# bin %% 2 means bin mod 2
      bin <- bin %/% 2 # this means bin integer division 2
      r})
    tune.mat <- replace(tune.mat.temp, !tune.mat.temp, l.forget) #replaces the values in tune.mat.temp with indices !tune.mat.temp with
                                                                #values of l.forget
    tune.mat <- replace(tune.mat, tune.mat.temp, l.noforget) #creates 2^d by d matrix with alternating rows of 'no forgetting' then
                                                              #'some forgetting' then 'no forgetting' then 'some forgetting etc.
    return(tune.mat)
  }

#-----------------------------------------------------------
## Updated: dlogr.step

dlogr.step <-
  function (x.t, y.t, betahat.tm1, varbetahat.tm1, tune.mat) 
  {
    if (!is.matrix(x.t)) {
      dim(x.t) <- c(1, length(x.t))
    }
    temp <- apply(tune.mat, 1, laplace.fn, x.t = x.t, y.t = y.t,  ##apllies the function laplace.fn to all rows of matrix tune.mat
                  betahat.tm1 = betahat.tm1, varbetahat.tm1 = varbetahat.tm1) #with the remaining variables passed to the function laplace.fn
    lambda <- tune.mat[1, ] #choose the value of lambda which maximises equation (8) in McCormick et al (2012)
    Rhat.t <- varbetahat.tm1
    Rhat.t <- (Rhat.t + t(Rhat.t))/2
    #print(lambda)
    diag(Rhat.t) <- diag(Rhat.t)/lambda
    laplace.t = max(temp)
    yhat.t <- dlogr.predict(x.t, betahat.tm1) #prediction of y at time t using x.t and betahat at time 1
    Del1 <- t(x.t) %*% (y.t - yhat.t)
    Del2 <- -solve(Rhat.t) - (t(x.t) * matrix(rep(yhat.t * (1 - yhat.t), dim(x.t)[2]), nrow = dim(x.t)[2], byrow = TRUE)) %*% 
              x.t
    betahat.t <- betahat.tm1 - (solve(Del2) %*% Del1)
    varbetahat.t <- solve(-Del2)
    diag(varbetahat.t) <- abs(diag(varbetahat.t))
    return(list(betahat.t = betahat.t, varbetahat.t = varbetahat.t, 
                laplace.t = laplace.t, lambda = lambda))
  }

    
#----------------------------------------------------------
##laplace.fn function
laplace.fn <- 
function (tune.vec, x.t, y.t, betahat.tm1, varbetahat.tm1) 
{
  if (!is.matrix(x.t)) {
    dim(x.t) <- c(1, length(x.t))
  }
  
  Rhat.t <- varbetahat.tm1
  Rhat.t <- (Rhat.t + t(Rhat.t))/2
  diag(Rhat.t) <- diag(Rhat.t)/tune.vec
  yhat.t <- dlogr.predict(x.t, betahat.tm1)
  Del1 <- t(x.t) %*% (y.t - yhat.t)
  Del2 <- -solve(Rhat.t) - (t(x.t) * matrix(rep(yhat.t * (1 -yhat.t), dim(x.t)[2]), nrow = dim(x.t)[2], byrow = TRUE)) %*% 
    x.t

  #Del2 <- Del2 - 0.00001
  #print(round(Rhat.t),2)  

  betahat.t <- betahat.tm1 - (solve(Del2) %*% Del1)



  
  p.theta <- dmnorm(t(betahat.t), t(betahat.tm1), Rhat.t) #multivariate normal distribution
  
  p.y <- prod((exp(y.t * x.t %*% betahat.t))/(1 + exp(x.t %*% betahat.t))) #likelihood function of logistic regression
  
  return((((dim(x.t)[2])/2) * log(2 * pi)) + (0.5 * log(abs(det(ginv((1 * Del2)))))) + log(p.theta) + log(p.y))
  #ginv Calculates the Moore-Penrose generalized inverse of a matrix X.
}

#----------------------------------------------------
#dlogr.predict function
dlogr.predict <- 
  function (x, BetaHat) 
  {
    if (!is.matrix(x)) {
      dim(x) <- c(1, length(x))
    }
    return(1/(1 + exp(-x %*% BetaHat))) #return inverse logit- i.e. the predicted mortality
  }

#-------------------------------------------------------
#dlogr.init function
dlogr.init <- 
function (x, y) 
{
  MyData <- data.frame(x, y)
  MyModel <- glm(y ~ x -1, data = MyData, family = binomial(link = logit)) #remove intercept since x already includes this in matrix
  rm(MyData)
  return(list(BetaHat = coefficients(MyModel), VarBetaHat = vcov(MyModel))) #returns estimated coefficients and var/cov matrix
}