library(RTMB)
y = readRDS("ML_ex/ar1.Rds")
plot(y,type = 'l',lwd = 2,cex.lab = 1.6, cex.axis = 1.5)

#Set up and estimate model
dat = list(y = y)
par = list(logsd = 0,
           logitRho = 0)

f<-function(par){
  getAll(dat,par)
  rho = 2/(1 + exp(-logitRho))-1
  timeSteps <-length(y)
  sd <- exp(logsd)
  
  nll<-0

  nll <- nll -dnorm(y[1],0,sqrt(sd*sd/(1-rho*rho)),log=TRUE)
  for(i in 2:timeSteps){    
    nll <- nll -dnorm(y[i],rho*y[i-1],sd,log=TRUE)
  }

  ADREPORT(rho) 
  nll
}

obj = MakeADFun(f,par)
opt = nlminb(obj$par,obj$fn,obj$gr)

AIC = 2*length(par) + 2*opt$objective

sdrep = sdreport(obj)
pl = as.list(sdrep,what = "Est")
plsd = as.list(sdrep,what = "Std")

# 95%CI interval for rho
c(2/(1 + exp(-pl$logitRho))-1,
  2/(1 + exp(-pl$logitRho + 2*plsd$logitRho*c(1,-1)))-1)


#Remove parameter with MAP-functionality



