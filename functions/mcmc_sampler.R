##
##
##
##
## functions for MCMC for all three models
## --------------------------------------

transformToGEV <- function(y,m,Qind,beta,xi,Delta){
  # transforms the index flood parameters to 
  # duration-dependent GEV parameters.
  # ----
  # outputs a 3 x n matrix of parameters where each column is 
  # a set of GEV params for a single duration
  # ---------------
  # m = number of years of data
  # y = data matrix with dim. m*n (num yrs * num durations) x 2
  # Qind,beta,xi = scalars
  # Delta = either a 1 x 1 (scalar) or a 1 x 2 vector of params
  #         depending on which model is running
  # ----------------------------------------
  n <- length(unique(y[,2])) #number of durations
  gevp <- matrix(ncol = n*3,nrow = 2) #GEV has 3 parameters
  colnames(gevp) <- rep(c("mu","sigma","xi"),n)
  j=1
  for(i in seq(1,m*n,by=m)){ #iterate over duration
    subY <- y[i:(i+(m-1)),]
    d <- unique(subY[,2])
    ## -- create the GEV parameters ----
    QindStar = Qind/(1+d*Delta[1])
    if(length(Delta)==2){ #check if we are in Javelle or DD
      sigma = QindStar*exp(beta)/(1+d*Delta[2])
    } else{sigma = QindStar*exp(beta)}
    mu = ifelse(abs(xi)>1e-15, 
                QindStar - sigma*(log(2)^(-xi)-1)/xi, 
                QindStar + sigma*log(log(2)))
    ## -- store the parameters ----
    gevp[,j:(j+2)] <- rbind(c(mu, sigma, xi),rep(d,3))
    j=j+3
  }
  return(gevp)
}

checkConstraints <- function(y,m,gevp,Delta){
  # checks constraints on:
  # (1) the Delta(s)
  # (2) the 3 x n GEV parameters 
  # ---------------
  # m = number of years of data
  # y = data matrix with dim. m*n (num yrs * num durations) x 2
  # gevp = matrix of gev parameters with dim. 3 x n 
  # Delta = either a 1 x 1 or a 1 x 2 vector of Delta params
  #         depending on which model is running
  # ----------------------------------------
  check = TRUE
  n <- length(unique(gevp[2,])) #number of durations
  yInd <- c(seq(1,m*n,by=m),(m*n+1)) #index to split data vector by duration
  pInd <- c(seq(1,n*3,by=3),(n*3+1)) #index to split parameter vec by duration
  ## -- Delta checks ----
  if(Delta[1]<=0){check = FALSE} #Deltas cannot be < 0
  if(length(Delta)==2){
    if(Delta[1]<Delta[2]){check = FALSE} #Delta2 cannot be > Delta1
    if(Delta[2]<=0){check = FALSE}} #Deltas cannot be < 0
  ## -- GEV support checks ----
  for(i in 1:n){
    subP <- gevp[1,pInd[i]:(pInd[i+1]-1)] #subset parameters by duration
    subY <- y[yInd[i]:(yInd[i+1]-1),1] #subset data by duration
    mu=subP["mu"]; sigma=subP["sigma"]; xi=subP["xi"]
    if(sigma <= 0){check = FALSE; break}
    if( any((1 + xi*(subY-mu)/sigma) <= 0) ){check = FALSE; break}
  }
  return(check)
}

loglik <- function(y,m,Qind,beta,xi,Delta){
  # 1. transform to GEV parameters
  gevp = transformToGEV(y,m,Qind,beta,xi,Delta) 
  # 2. check constraints
  if(checkConstraints(y,m,gevp,Delta) != TRUE){return(-Inf)}
  # 3. sum the log likelihood over both durations and data points
  total = 0
  n <- length(unique(gevp[2,]))
  yInd <- c(seq(1,m*n,by=m),(m*n+1)); pInd <- c(seq(1,n*3,by=3),(n*3+1)) 
  for(i in 1:n){ #loop over likelihood for each duration
    subP <- gevp[1,pInd[i]:(pInd[i+1]-1)]; subY <- y[yInd[i]:(yInd[i+1]-1),1]
    mu=subP["mu"]; sigma=subP["sigma"]; xi=subP["xi"]
    total = total + sum(extRemes::devd(subY,mu,sigma,xi,
                                       log = T,type = 'GEV'))}
  return(total)
}

logprior <- function(y,m,Qind,beta,xi,Delta){
  # y, m passed in solely for constraint checking purposes
  gevp = transformToGEV(y,m,Qind,beta,xi,Delta) 
  if(checkConstraints(y,m,gevp,Delta) != TRUE){return(-Inf)}
  # construct log prior:
  Qpr = log(truncnorm::dtruncnorm(Qind,a=0,b=Inf,40,100))
  bpr = dnorm(beta, 0, 100, log = T)
  xpr = dbeta(-xi + 0.5, shape1=6, shape2=9, log=TRUE)
  if(length(Delta)==2){
    Dpr = dlnorm(Delta[2],0,5,log = T)+log(dlnormTrunc(Delta[1],0,5,Delta[2]))
  }else{Dpr = dlnorm(Delta[1], 0, 5, log = T)}
  return(Qpr+bpr+xpr+Dpr)
}

checkStartpoints <- function(y,tuning,m,Qind,beta,xi,Delta){
  gevp = transformToGEV(y,m,Qind,beta,xi,Delta) 
  if(checkConstraints(y,m,gevp,Delta) != TRUE){
    check = FALSE # if false, propose new startpoints
    while(check != TRUE){
      Qnew=rnorm(1,mean=Qind,sd=2*tuning["qIndSD"])
      Bnew=runif(1,beta-0.15,beta+0.15)
      Xnew=runif(1,xi-0.15,xi+0.15)
      if(length(Delta)==2){
        Dnew=c(rnorm(1,mean=Delta[1],sd=2*tuning["Delta1SD"]),
               rnorm(1,mean=Delta[2],sd=2*tuning["Delta2SD"]))
      }else{Dnew=rnorm(1,mean=Delta,sd=2*tuning["DeltaJSD"])}
      gevp = transformToGEV(y,m,Qnew,Bnew,Xnew,Dnew) 
      check = checkConstraints(y,m,gevp,Dnew)}
    Qind=Qnew; beta=Bnew; xi=Xnew; Delta=Dnew
  }
  return(list(Qind,beta,xi,Delta))
}

javelle <- function(data,startpoint,tuning,iter,ss){
  m = dim(data)[1]/length(unique(data[,2]))
  stpt=checkStartpoints(data,tuning,m,startpoint["Q"],startpoint["B"],
                        startpoint["X"],startpoint["D1"])
  Qind=stpt[[1]];beta=stpt[[2]];xi=stpt[[3]];Delta=stpt[[4]]
  oll=loglik(data,m,Qind,beta,xi,Delta)
  olpr=logprior(data,m,Qind,beta,xi,Delta)
  
  Nsamp = iter/ss; 
  obj=matrix(NA,Nsamp,3); colnames(obj)<-c("numDeltas","logprior","loglik")
  accept=rep(0,4); names(accept) <- c("Qind","beta","xi","Delta")
  params=matrix(NA,Nsamp,4); colnames(params) <- c("Qind","beta","xi",
                                                   "Delta")
  for(i in 1:iter){
    # -- update Qind -----
    Qnew=rnorm(1,mean=Qind,sd=2*tuning["qIndSD"])
    nll=loglik(data,m,Qnew,beta,xi,Delta)
    nlpr=logprior(data,m,Qnew,beta,xi,Delta)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qnew; beta=beta; xi=xi; Delta=Delta
      oll=nll; olpr=nlpr
      accept[1]=accept[1]+1}
    # --------------------
    # -- update beta -----
    Bnew=runif(1,beta-0.05,beta+0.05)
    nll=loglik(data,m,Qind,Bnew,xi,Delta)
    nlpr=logprior(data,m,Qind,Bnew,xi,Delta)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qind; beta=Bnew; xi=xi; Delta=Delta
      oll=nll; olpr=nlpr
      accept[2]=accept[2]+1}
    # --------------------
    # -- update xi -----
    Xnew=runif(1,xi-0.05,xi+0.05)
    nll=loglik(data,m,Qind,beta,Xnew,Delta)
    nlpr=logprior(data,m,Qind,beta,Xnew,Delta)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qind; beta=beta; xi=Xnew; Delta=Delta
      oll=nll; olpr=nlpr
      accept[3]=accept[3]+1}
    # --------------------
    # -- update Delta -----
    Dnew=rnorm(1,mean=Delta,sd=2*tuning["DeltaJSD"])
    nll=loglik(data,m,Qind,beta,xi,Dnew)
    nlpr=logprior(data,m,Qind,beta,xi,Dnew)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qind; beta=beta; xi=xi; Delta=Dnew
      oll=nll; olpr=nlpr
      accept[4]=accept[4]+1}
    # --------------------
    # store every ss-th iter
    if (i%%ss==0) {
      obj[i/ss,]=c(length(Delta),olpr,oll)
      params[i/ss,]=c(Qind,beta,xi,Delta)}
  }
  return(list(obj,params,accept))
}


extendedQDF <- function(data,startpoint,tuning,iter,ss){
  m = dim(data)[1]/length(unique(data[,2]))
  stpt=checkStartpoints(data,tuning,m,startpoint["Q"],
                        startpoint["B"],startpoint["X"],
                        c(startpoint["D1"],startpoint["D2"]))
  Qind=stpt[[1]];beta=stpt[[2]];xi=stpt[[3]];Delta=stpt[[4]]
  oll=loglik(data,m,Qind,beta,xi,Delta)
  olpr=logprior(data,m,Qind,beta,xi,Delta)
  
  Nsamp = iter/ss; 
  obj=matrix(NA,Nsamp,3); colnames(obj)<-c("numDeltas","logprior","loglik")
  accept=rep(0,5); names(accept) <- c("Qind","beta","xi","Delta1","Delta2")
  params=matrix(NA,Nsamp,5); colnames(params) <- c("Qind","beta","xi",
                                                   "Delta1","Delta2")
  for(i in 1:iter){
    # -- update Qind -----
    Qnew=rnorm(1,mean=Qind,sd=2*tuning["qIndSD"])
    nll=loglik(data,m,Qnew,beta,xi,Delta)
    nlpr=logprior(data,m,Qnew,beta,xi,Delta)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qnew; beta=beta; xi=xi; Delta=Delta
      oll=nll; olpr=nlpr
      accept[1]=accept[1]+1}
    # --------------------
    # -- update beta -----
    Bnew=runif(1,beta-0.05,beta+0.05)
    nll=loglik(data,m,Qind,Bnew,xi,Delta)
    nlpr=logprior(data,m,Qind,Bnew,xi,Delta)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qind; beta=Bnew; xi=xi; Delta=Delta
      oll=nll; olpr=nlpr
      accept[2]=accept[2]+1}
    # --------------------
    # -- update xi -----
    Xnew=runif(1,xi-0.05,xi+0.05)
    nll=loglik(data,m,Qind,beta,Xnew,Delta)
    nlpr=logprior(data,m,Qind,beta,Xnew,Delta)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qind; beta=beta; xi=Xnew; Delta=Delta
      oll=nll; olpr=nlpr
      accept[3]=accept[3]+1}
    # --------------------
    # -- update Deltas -----
    Dnew=c(rnorm(1,mean=Delta[1],sd=2*tuning["Delta1SD"]),
           Delta[2])
    nll=loglik(data,m,Qind,beta,xi,Dnew)
    nlpr=logprior(data,m,Qind,beta,xi,Dnew)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qind; beta=beta; xi=xi; Delta=Dnew
      oll=nll; olpr=nlpr
      accept[4]=accept[4]+1}

    Dnew=c(Delta[1],
           rnorm(1,mean=Delta[2],sd=2*tuning["Delta2SD"]))
    nll=loglik(data,m,Qind,beta,xi,Dnew)
    nlpr=logprior(data,m,Qind,beta,xi,Dnew)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qind; beta=beta; xi=xi; Delta=Dnew
      oll=nll; olpr=nlpr
      accept[5]=accept[5]+1}
    # --------------------
    # store every ss-th iter
    if (i%%ss==0) {
      obj[i/ss,]=c(length(Delta),olpr,oll)
      params[i/ss,]=c(Qind,beta,xi,Delta)}
  }
  return(list(obj,params,accept))
}

reversiblejumpQDF <- function(data,startpoint,tuning,iter,ss,innerloop){
  m = dim(data)[1]/length(unique(data[,2]))
  stpt=checkStartpoints(data,tuning,m,startpoint["Q"],
                        startpoint["B"],startpoint["X"],
                        c(startpoint["D1"],startpoint["D2"]))
  Qind=stpt[[1]];beta=stpt[[2]];xi=stpt[[3]];Delta=stpt[[4]]
  oll=loglik(data,m,Qind,beta,xi,Delta)
  olpr=logprior(data,m,Qind,beta,xi,Delta)
  
  Nsamp = iter/ss; 
  obj=matrix(NA,Nsamp,3); colnames(obj)<-c("numDeltas","logprior","loglik")
  params=list()
  for(i in 1:iter){
    for(j in 1:innerloop){
    # -- update Qind -----
    Qnew=rnorm(1,mean=Qind,sd=2*tuning["qIndSD"])
    nll=loglik(data,m,Qnew,beta,xi,Delta)
    nlpr=logprior(data,m,Qnew,beta,xi,Delta)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qnew; beta=beta; xi=xi; Delta=Delta
      oll=nll; olpr=nlpr}
    # --------------------
    # -- update beta -----
    Bnew=runif(1,beta-0.05,beta+0.05)
    nll=loglik(data,m,Qind,Bnew,xi,Delta)
    nlpr=logprior(data,m,Qind,Bnew,xi,Delta)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qind; beta=Bnew; xi=xi; Delta=Delta
      oll=nll; olpr=nlpr}
    # --------------------
    # -- update xi -----
    Xnew=runif(1,xi-0.05,xi+0.05)
    nll=loglik(data,m,Qind,beta,Xnew,Delta)
    nlpr=logprior(data,m,Qind,beta,Xnew,Delta)  
    ratio=nll+nlpr-oll-olpr
    if (log(runif(1))<ratio) {
      Qind=Qind; beta=beta; xi=Xnew; Delta=Delta
      oll=nll; olpr=nlpr}
    # --------------------
    # -- update Deltas -----
    if(length(Delta)==2){ #update following double-delta
      Dnew=c(rnorm(1,mean=Delta[1],sd=2*tuning["Delta1SD"]),
             Delta[2])
      nll=loglik(data,m,Qind,beta,xi,Dnew)
      nlpr=logprior(data,m,Qind,beta,xi,Dnew)  
      ratio=nll+nlpr-oll-olpr
      if (log(runif(1))<ratio) {
        Qind=Qind; beta=beta; xi=xi; Delta=Dnew
        oll=nll; olpr=nlpr}
      Dnew=c(Delta[1],
             rnorm(1,mean=Delta[2],sd=2*tuning["Delta2SD"]))
      nll=loglik(data,m,Qind,beta,xi,Dnew)
      nlpr=logprior(data,m,Qind,beta,xi,Dnew)  
      ratio=nll+nlpr-oll-olpr
      if (log(runif(1))<ratio) {
        Qind=Qind; beta=beta; xi=xi; Delta=Dnew
        oll=nll; olpr=nlpr}
    }else{ #update following Javelle
      Dnew=rnorm(1,mean=Delta,sd=2*tuning["DeltaJSD"])
      nll=loglik(data,m,Qind,beta,xi,Dnew)
      nlpr=logprior(data,m,Qind,beta,xi,Dnew)  
      ratio=nll+nlpr-oll-olpr
      if (log(runif(1))<ratio) {
        Qind=Qind; beta=beta; xi=xi; Delta=Dnew
        oll=nll; olpr=nlpr}}
    # -----------------------------
    } #end of inner loop
    # -- split / combine step -----
    if(length(Delta)==1){ #try to split Delta
      u <- rbeta(1,5,1)
      splitD <- c(Delta*u,(1-u)*Delta)
      q <- dbeta(u,5,1,log=T)
      J <- log(Delta)
      nll=loglik(data,m,Qind,beta,xi,splitD)
      nlpr=logprior(data,m,Qind,beta,xi,splitD)
      ratio=nll+nlpr-oll-olpr-q+J
      if (log(runif(1))<ratio) {
        Qind=Qind; beta=beta; xi=xi; Delta=splitD
        oll=nll; olpr=nlpr}
    }else{ #try to combine Delta
      combD <- Delta[1] + Delta[2]
      q <- dbeta(Delta[1]/(Delta[1]+Delta[2]),5,1,log=T)
      J <- log(1/(Delta[1]+Delta[2]))
      nll=loglik(data,m,Qind,beta,xi,combD)
      nlpr=logprior(data,m,Qind,beta,xi,combD)
      ratio=nll+nlpr-oll-olpr+q+J
      if (log(runif(1))<ratio) {
        Qind=Qind; beta=beta; xi=xi; Delta=combD
        oll=nll; olpr=nlpr}}
    # --------------------
    # store every ss-th iter
    if (i%%ss==0) {
      obj[i/ss,]=c(length(Delta),olpr,oll)
      params[[i/ss]]=list(Qind=Qind,beta=beta,xi=xi,Delta=Delta)}
  }
  return(list(obj,params))
}