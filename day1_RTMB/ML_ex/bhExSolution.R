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

  sigma = exp(logSigma)
  
  pred_logR = logA + log(SSB) -log(1 + exp(logB)*SSB)
  nll = -sum(dnorm(logR, pred_logR,sigma,TRUE))

  ADREPORT(pred_logR)  
  return(nll);
}


obj <- MakeADFun(f,par)
opt <- nlminb(obj$par,obj$fn,obj$gr)
sdrep = sdreport(obj)

rl = as.list(sdrep,what = "Est",report = TRUE)
rlsd = as.list(sdrep,what = "Std",report = TRUE)

plot(dat$SSB,dat$logR)
lines(dat$SSB, rl$pred_logR)
lines(dat$SSB,rl$pred_logR + 2*rlsd$pred_logR,lty = 2)
lines(dat$SSB,rl$pred_logR - 2*rlsd$pred_logR,lty = 2)


