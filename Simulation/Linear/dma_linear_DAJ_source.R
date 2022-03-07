library(dma)

dma <-
  function(x, y, models.which, lambda=0.99, gamma=0.99, 
           eps=.001/nrow(models.which), delay=0, initialperiod=200, update_t) UseMethod("dma")

dma.default <-
  function(x, y, models.which, lambda=0.99, gamma=0.99, 
           eps=.001/nrow(models.which), delay=0, initialperiod=200, update_t) {
    
    #run 
    est<-makf4(x, y, models.which, lambda, gamma,eps, delay, initialperiod, update_t)
    
    #define outputs
    est$fitted.values<-est$yhat.ma
    est$residuals<-y-est$yhat.ma
    
    #define as class to allow other commands later
    class(est)<-"dma"
    
    est
  }

rm.Kalman <-
  function (thetaold, Sigmaold, Vold, lambda, x, y, tt) {
    
    # August 8, 2007. Function to do Kalman prediction and updating
    #  for a specific model in the rolling mill problem.
    #  Here V_t==V, G_t==1,  and the state equation is given by forgetting.
    
    # Inputs:
    # thetaold	vector containing initial state estimate, theta hat_{t-1}
    # Sigmaold	square matrix containing initial state covariance, Sigma_{t-1}
    # Vold		innovations variance
    # lambda	forgetting factor for state equation
    # x		x_t = system input at time t
    # y		y_t = system output at time t (1-dimensional)
    # tt		index of new observation (t) - needed to update V
    
    # Outputs:
    # thetanew	updated state estimate, theta hat_t
    # Sigmanew	posterior covariance matrix of updated state estimate
    # yhat		Predictive mean of y_t: E[y_t|x,y^{t-1}]
    # predvar	Predictive variance of y_t: Var[y_t|x,y^{t-1}]
    # Vnew		new estimate of innovations variance V
    
    # Initial calculations
    x <- matrix(x, ncol=dim(x)[2])
    #print(x)
    Vold1 <- Vold
    # Kalman updating
    R <- Sigmaold/lambda					# Rt

    Vold <- diag(Vold,length(x[,1]))
    #K <- R %*% x %*%  solve( Vold + t(x) %*% R %*% x)	# Kalman gain
    K <- (R %*% t(x)) %*%  solve( Vold + (x %*% R %*% t(x)))	# Kalman gain matrix
    #K <- (R %*% t(x)) %*%  solve( Vold1 + mean(diag((x %*% R %*% t(x)))))	# Kalman gain matrix
    yhat <- x %*% thetaold				# Predictive mean
    predvar <- (x %*% R %*% t(x)) + Vold			# Predictive variance
    #predvar <- diag(c(diag(x %*% R %*% t(x))),dim(x)[1]) + Vold			# Predictive variance
    e <- y - yhat 						# Prediction error
    thetanew <- thetaold + (K %*% e)
    #print(K %*% e)
    #print(K)
    #print(e)
    #print(e)
    Sigmanew <- R - K %*% x %*% R
    
    # Updating V; see C8303
    #### CSP added -1 in line below to get same estimates as original code, changed 'sum' to 'mean'
    Vnew <- ((min(tt)-1)/max(tt)) * Vold1 + ((max(tt)-min(tt)+1)/max(tt))*mean(diag(e %*% t(e) - x %*% R %*% t(x)))
    td1 <- e %*% t(e) - x %*% R %*% t(x)
    #print(paste0("V old=",Vold1,", V new=",Vnew))
    #print(paste0("theta old=",thetaold,", theta new=",thetanew))
    if (Vnew <= 0) Vnew <- Vold1
    
    # Output
    list (thetanew=thetanew, Sigmanew=Sigmanew, yhat=yhat, predvar=mean(diag(predvar)), 
          Vnew=Vnew)
  }

makf4 <- function( x, y, models.which, lambda=0.99, gamma=0.99, 
                   eps=.001/nrow(models.which), delay=0, initialperiod=200, update_t) {
  # Revised June 29, 2009 to replace equation (17) in the paper by a version
  # (Minor revision of comments August 10, 2011.)
  #  that regularizes all posterior model probabilities away from zero.
  #  This is to avoid degenercies in which models are eliminated for ever,
  #  or model probabilities go to machine zero.
  # Revised June 25, 2009 to output quantities needed for simulations.
  # August 23, 2007.
  # Function to implement the Model Averaging Kalman Filter (MAKF) for
  #  the rolling mill problem. See C8300-C8303 and C3883.
  # In this version, the model probabilities are updated by flattening.
  # Calls functions rm.Kalman and model.update3.
  # This version is a revised version of makf, revised to take account of the
  #  fact that when yt is being controlled, only information on y_{t-delay-1} and 
  #  earlier is available for estimating theta. See C8317-8.
  
  # Inputs:
  # x		TTxp matrix of system inputs
  # y		TT-vector of system outputs
  #  models.which	Kxp matrix, with 1 row per model and 1 col per variable
  #		indicating whether that variable is in the model
  #  (the state theta is of dim (model.dim+1); the extra 1 for the intercept)
  # lambda	parameter forgetting factor
  # gamma		flatterning parameter for model updating
  # eps		regularization parameter for regularizing posterior model
  #		 model probabilities away from zero
  # delay		When y_t is controlled, only y_{t-delay-1} and before 
  #		 are available. This is determined by the machine.
  #		Note that delay as defined here corresponds to (k-1)
  #		 in the Ettler et al (2007, MixSim) paper.
  #		Thus we use the default delay=24, corresponding to k=25.
  # initialperiod length of initial period. Performance is summarized with and
  #		 without the first initialperiod samples.
  
  # Outputs
  # yhat.bymodel	TTxK matrix whose tk element gives yhat for yt for model k
  # yhat.ma	TT vector whose t element gives the model-averaged yhat for yt
  # pmp		TTxK matrix whose tk element is the post prob of model k at t
  # thetahat.ma	TTx(nvar+1) matrix whose tk element is the model-averaged 
  #		 estimate of theta_{j-1} at t
  # Vtheta.ma	TTx(nvar+1) matrix whose tk element is the model-averaged
  #		 variance of thetahat_{j-1} at t
  # mse.bymodel	MSE for each model
  # mse.ma	MSE of model-averaged prediction
  # mse400.bymodel  MSE for each model excluding the first initialperiod samples
  # mse400.ma	MSE of model averaging excluding the first initialperiod samples
  # model.forget	forgetting factor for the model switching matrix
  
  # Initialization
  x <- as.matrix (x)
  nvar <- ncol(x)
  TT <- length(y)
  K <- nrow (models.which)
  model.dim <- rep(0,K)
  for (k in 1:K) model.dim[k] <- sum(as.numeric(models.which[k,]))+1
  thetaold <- vector ("list",K)
  Sigmaold <- vector ("list",K)
  yhat <- matrix (rep(0,TT*K), ncol=K)
  if (delay>0) yhat[1:delay,] <- NA
  yhat.ma <- rep( 0,TT)
  if (delay>0) yhat.ma[1:delay] <- NA
  yhat.Kalman <- matrix (rep(0,TT*K), ncol=K)
  pimat <- matrix (rep(0,(TT-delay)*K), ncol=K)
  
  # Added on 6/25/09
  thetahat <- matrix (rep(0,K*(nvar+1)),nrow=K)  # temporary variable
  Vtheta <- matrix (rep(0,K*(nvar+1)),nrow=K)  # temporary variable
  thetahatall <- Vthetaall <- array(0, dim=c(K,TT,nvar+1))
  thetahat.ma <- matrix (rep(0,TT*(nvar+1)),nrow=TT)
  Vtheta.ma <- matrix (rep(0,TT*(nvar+1)),nrow=TT)
  # End 6/25/09 addition
  
  # Initialize thetaold, Sigmaold, Vold, piold, for each model
  s2y <- 55.6 # variance of y
  Vold <- rep( s2y, K)
  s2x <- rep(0,nvar)
  for (j in (1:nvar)) s2x[j] <- var(x[,j])
  for (k in 1:K) {
    thetaold[[k]] <- c (0, rep( 0, model.dim[k]-1) )
    if (model.dim[k] == 1) Sigmaold[[k]] <- 430^2 else
      Sigmaold[[k]] <- diag( c( 430^2 , s2y/s2x[as.logical(models.which[k,])]) )
    # The intercept prior mean and variance are UIP for the whole dataset
  }
  piold <- rep (1/K, K)
  
  # Iterations
  
  for (tx in 2:length(update_t)) {
    predvar <- rep(0,K)
    
    t1 <- update_t[tx-1]+1
    t2 <- update_t[tx]
    tt <- t1:t2
    
    # Kalman update for each model
    for (k in 1:K) {
      which <- as.logical( models.which[k,])
      xtmp1 <-  matrix(x[tt,which], ncol=dim(x)[2])
      
      xtmp <-  cbind(1,xtmp1)
      
      Kalman.out <- rm.Kalman (thetaold[[k]], Sigmaold[[k]], Vold[k], lambda,
                               xtmp, y[tt], tt)
      
      yhat.Kalman[tt,k] <- Kalman.out$yhat
      #ta <- matrix(thetaold[[k]],nrow=(max(tt)-min(tt)),ncol=length(thetaold[[k]]),byrow=TRUE)
      #print(dim(ta))
      #print(dim(t(c(1,x[tt+delay,which])))) # needs changing to matric
      #print(dim(t(c(x[tt+delay,which]))))
      #yhat[tt+delay,k] <- t(c(1,x[tt+delay,which])) %*% ta
      # Update theta and Sigma using Kalman filter output
      thetaold[[k]] <- Kalman.out$thetanew
      Sigmaold[[k]] <- Kalman.out$Sigmanew
      #print(Kalman.out$Vnew)
      Vold[k] <- Kalman.out$Vnew
      predvar[k] <- Kalman.out$predvar
      
      # Added on 6/25/09
      whichtmp <- c(TRUE,which)
      thetahat[k,whichtmp] <- Kalman.out$thetanew
      Vtheta[k,whichtmp] <- diag (Kalman.out$Sigmanew)
      # End 6/25/09 addition
    }
    
    # Model averaging prediction
    yhat.ma[tt+delay] <- t(piold) %*% yhat[tt+delay,]
    
    # Update model probabilities
    model.update.out <- model.update3 (piold, gamma, eps, y[tt], yhat.Kalman[tt,], predvar)
    piold <- model.update.out$pinew
    pimat[tt,] <- piold
    
    # Added 6/25/09: Model averaged estimates and posterior variances
    thetahat.ma[tt,] <- t(piold) %*% thetahat
    
    
    
    #Vtheta.ma[tt,] <- t(piold) %*% Vtheta + t(piold) %*% (thetahat %*% t(thetahat)) - (thetahat.ma[tt,] %*% t(thetahat.ma[tt,]))
    
    #th1 <- matrix(thetahat,nrow=(max(tt)-min(tt)),ncol=length(thetahat),byrow=TRUE)
    #print(thetahat)
    #thetahatall[, tt, ] <- th1
    for(j in min(tt):max(tt)){
      thetahatall[, j, ] <- thetahat
    }
    
    Vthetaall[, tt, ] <- Vtheta
  }
  
  # Summarize performance with and without the first initialperiod samples.
  #  (If initialperiod < delay+1, then initialperiod <- (delay+1) is used.
  #   In that case, performance with and without the first initialperiod
  #   samples should be the same.)
  # The first observation (or the (delay+1)-th observation) is ignored
  #  because its prediction is always zero.
  # The measures without the first initialperiod samples are called 
  #  mse400.bymodel and mse400.ma, because Petr Nedoma used initialperiod=400.
  initialperiod <- max(delay+1,initialperiod)
  mse.bymodel <- rep(0,K)
  mse400.bymodel <- rep(0,K)
  for (k in 1:K) mse.bymodel[k] <- mean ( (y[(delay+2):TT] - yhat[(delay+2):TT,k])^2 ) 
  for (k in 1:K) mse400.bymodel[k] <- 
    mean ( (y[(initialperiod+1):TT] - yhat[(initialperiod+1):TT,k])^2 )
  mse.ma <- mean ( (y[(delay+2):TT]-yhat.ma[(delay+2):TT])^2 )
  mse400.ma <- mean ( (y[(initialperiod+1):TT] - yhat.ma[(initialperiod+1):TT])^2 )
  msemseinitialperiod.bymodel<-mse400.bymodel
  msemseinitialperiod.ma<-mse400.ma
  
  # Output
  list (yhat.bymodel=yhat, yhat.ma=yhat.ma, pmp=pimat, thetahat=thetahatall, Vtheta=Vthetaall,
        thetahat.ma=thetahat.ma, Vtheta.ma=Vtheta.ma,
        mse.bymodel=mse.bymodel, mse.ma=mse.ma, mseinitialperiod.bymodel=msemseinitialperiod.bymodel,
        mseinitialperiod.ma=msemseinitialperiod.ma, lambda=lambda, gamma=gamma,
        models.which=models.which, delay=delay, initialperiod=initialperiod)
}

model.update3 <-
  function (piold, gamma, eps,  y, yhat, predvar) {
    
    # Revised June 29, 2009:
    # Modified to regularize the posterior model probabilities away from zero
    #  by adding eps to each one and renormalizing.
    # August 23, 2007. Update model posterior probabilities using 
    #  flattening. See C8338-9.
    # This will be used in makf3.
    
    # Inputs:
    #  piold	K-vector of input model probabilities
    #  gamma	flattening parameter
    #  eps		minimum threshold for model probabilities  
    #  y		observed value of y_t
    #  yhat		K-vector of predicted values of y_t | y_{t-1} from rm.Kalman
    #  predvar	K-vector of predicted variances of y_t | y_{t-1} from rm.Kalman
    
    # Output:
    #  pinew	K-vector of updated model probabilities
    
    # Form predicted pi values
    pipred <- piold^gamma / sum(piold^gamma)
    
    # Update pi values
    logpyt <- -0.5*log(predvar) - 0.5*(y-yhat)^2/predvar
    logpyt <- logpyt - max(logpyt)
    pyt <- exp (logpyt)
    pinew <- pipred * pyt
    pinew <- pinew/sum(pinew)
    pinew <- pinew + eps
    pinew <- pinew/sum(pinew)
    pinew <- mean(pinew)
    
    # Output
    list (pinew=as.vector(pinew))
  }

