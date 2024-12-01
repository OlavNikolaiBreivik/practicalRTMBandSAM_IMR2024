load("introStockAssessment_ex/fsa.RData") # gets "dat"

library(RTMB)

dat$keyFsta =  c(1,2,3,4,5,5,5)
par <- list(
  logN = dat$M*0,
  logF = dat$M[1:max(dat$keyFsta),],
  sdF = 0,
  sdN = c(0,0),
  logSdCatch=0,
  logQ=rep(0,length(unique(dat$age[dat$fleet==2]))),
  logSdSurvey=0
)

f<-function(par){
  getAll(par, dat)
  nll = 0

  na <- max(age)-min(age)+1
  ny <- max(year)-min(year)+1

  #Fishing mortality
  SigmaF = diag(dim(logF)[1]) * exp(sdF)^2
  for(i in 2:dim(logF)[2]){
    nll = nll- dmvnorm(logF[,i], logF[,i-1],SigmaF, TRUE) #Random walk for logF
  }
  ## setup F
  F <- exp(logF[keyFsta,])

  ## setup N
  for(y in 2:ny){
    nll = nll - dnorm(logN[1,y],logN[1,y-1] , exp(sdN[1]), TRUE) #Random walk for recruitment
    for(a in 2:na){
      nll = nll- dnorm(logN[a,y],logN[a-1,y-1]-F[a-1,y-1]-M[a-1,y-1], exp(sdN[2]),TRUE) #Survival process
    }
  }

  # Match to observations
  logObs <- log(obs)
  logPred <- numeric(length(logObs))
  sdvec <- numeric(length(logObs))
  for(i in 1:length(logObs)){
    a <- age[i]-min(age)+1
    y <- year[i]-min(year)+1
    if(fleet[i]==1){
      logPred[i] <- log(F[a,y])-log(F[a,y]+M[a,y])+log(1-exp(-F[a,y]-M[a,y]))+logN[a,y]
      sdvec[i] <- exp(logSdCatch)
    }else{
      logPred[i] <- logQ[a]-(F[a,y]+M[a,y])*surveyTime+logN[a,y]
      sdvec[i] <- exp(logSdSurvey)
    }
  }

  nll <-nll -sum(dnorm(logObs,logPred,sdvec,TRUE))

  logssb <- log(apply(exp(logN)*stockMeanWeight*propMature,2,sum))
  ADREPORT(logssb)
  return(nll)
}


#Set up AD-machinery
obj <- MakeADFun(f, par, random = c("logF","logN"))

#Estimate model
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(iter.max=1000,eval.max=1000, trace = 1))

#Extract adreported variables
sdrep <- sdreport(obj)


rlStateSpace <- as.list(sdrep, "Est", report=TRUE)
rlsdStateSpace <- as.list(sdrep, "Std", report=TRUE)

yr<-sort(unique(dat$year))
plot(yr, exp(rlStateSpace$logssb), type="l", lwd=5, col="red", ylim=c(0,550000), xlab="Year", ylab="SSB")
lines(yr, exp(rlStateSpace$logssb-2*rlsdStateSpace$logssb), type="l", lwd=1, col="red")
lines(yr, exp(rlStateSpace$logssb+2*rlsdStateSpace$logssb), type="l", lwd=1, col="red")



#Compare with the parametric assessment,
lines(yr, exp(pl$logssb), type="l", lwd=5, col="blue")
lines(yr, exp(pl$logssb-2*plsd$logssb), type="l", lwd=1, col="blue")
lines(yr, exp(pl$logssb+2*plsd$logssb), type="l", lwd=1, col="blue")


