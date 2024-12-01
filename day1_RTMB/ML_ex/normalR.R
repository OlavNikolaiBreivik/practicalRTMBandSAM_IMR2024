library(RTMB)

dat = list(Y = c(0.8,-1.0,-0.7, -2.0,-2.1,0.6,0.9,0.1, 1.8,4.1,-1.0))
par = list(mu = 0,
           logsd = 0)

f = function(par){
  getAll(dat,par)
  sd= exp(logsd)
  nll= -sum(dnorm(Y,mu,sd,log=TRUE))
  return(nll)
}

obj = MakeADFun(f,par)
opt = nlminb(obj$par,obj$fn,obj$gr, control = list(trace = 1))
