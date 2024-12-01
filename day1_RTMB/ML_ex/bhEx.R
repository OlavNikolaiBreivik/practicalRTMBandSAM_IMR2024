dat<-read.table("ML_ex/bh.dat", header=TRUE)
library(RTMB)


data <- list(SSB=dat$SSB,logR=dat$logR)
par <- list(
  logA=0,
  logB=0,
  logSigma=0
)


f = function(par){
  getAll(par,data)

  pred_logR = ...
  
  nll = -sum(dnorm(...))
  
  #Adreport predictions for plotting
  ADREPORT(pred_logR)
  
  return(nll);
}


obj <- MakeADFun(f,par)
opt <- nlminb(...)

