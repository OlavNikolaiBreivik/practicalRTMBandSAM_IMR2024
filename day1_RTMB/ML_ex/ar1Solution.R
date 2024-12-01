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

# 95%CI intervals
c(2/(1 + exp(-pl$logitRho))-1,
  2/(1 + exp(-pl$logitRho + 2*plsd$logitRho*c(1,-1)))-1)


#Remove parameter with MAP-functionality
par$logitRho = 0
map = list(logitRho = as.factor(NA))
obj2 = MakeADFun(f,par,map = map)
opt2 = nlminb(obj2$par,obj2$fn,obj2$gr)

AIC2 = 2*(length(par)-1) + 2*opt2$objective

#Liklihood-ratio test, more on this later.
LR = -2*(opt$objective -opt2$objective)
1-pchisq(LR,df = 1)



