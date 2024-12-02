library(stockassessment)

cn<-read.ices("nssherring/cn.dat")
cw<-read.ices("nssherring/cw.dat")
dw<-read.ices("nssherring/dw.dat")
lf<-read.ices("nssherring/lf.dat")
lw<-read.ices("nssherring/lw.dat")
mo<-read.ices("nssherring/mo.dat")
nm<-read.ices("nssherring/nm.dat")
pf<-read.ices("nssherring/pf.dat")
pm<-read.ices("nssherring/pm.dat")
sw<-read.ices("nssherring/sw.dat")
surveys<-read.ices("nssherring/survey.dat")

dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn, 
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

confDef<-defcon(dat)
par = defpar(dat,confDef)
fit = sam.fit(dat,confDef,par)

saveConf(confDef,file = "nssherring/conf.cfg")

fit$sdrep


res = residuals(fit)

residDiagPlot(fit,resid = res)

conf01Fsta = loadConf(dat,"nssherring/suggestedConf/conf01Fsta.cfg")
par01Fsta = defpar(dat,conf01Fsta)
fit01Fsta = sam.fit(dat,conf01Fsta,par01Fsta)
partable(fit01Fsta)


##########
conf02 = loadConf(dat,"nssherring/suggestedConf/conf02LogFpar.cfg")
par02 = defpar(dat,conf02)
fit02 = sam.fit(dat,conf02,par02)
partable(fit02LogFpar)

AIC(fit,fit01Fsta,fit02)
res02 = residuals(fit02)
residDiagPlot(fit02,resid = res02)

##########
conf03 = loadConf(dat,"nssherring/suggestedConf/conf03PredVar.cfg")
par03 = defpar(dat,conf03)
fit03 = sam.fit(dat,conf03,par03)
partable(fitTmp)

AIC(fit,fit01Fsta,fit02,fit03)
res03 = residuals(fit03)
residDiagPlot(fit03,resid = res03)


##########
conf04 = loadConf(dat,"nssherring/suggestedConf/conf04ObsVar.cfg")
par04 = defpar(dat,conf04)
fit04 = sam.fit(dat,conf04,par04)
partable(fit04)

AIC(fit,fit01Fsta,fit02,fit03,fit04)
res04 = residuals(fit04)
residDiagPlot(fitTmp,resid = res04)


##########
conf05 = loadConf(dat,"nssherring/suggestedConf/conf05CorObs.cfg")
par05 = defpar(dat,conf05)
fit05 = sam.fit(dat,conf05,par05)
partable(fit05)

AIC(fit,fit01Fsta,fit02,fit03,fit04,fit05)
resTmp = residuals(fitTmp)
residDiagPlot(fitTmp,resid = resTmp)


##########
confTmp = loadConf(dat,"nssherring/suggestedConf/confTmp.cfg")
parTmp = defpar(dat,confTmp)
fitTmp = sam.fit(dat,confTmp,parTmp)
partable(fitTmp)

AIC(fit,fit01Fsta,fit02,fit03,fit04,fit05,fitTmp)
resTmp = residuals(fitTmp)
residDiagPlot(fitTmp,resid = resTmp)




####
ssbplot(c(fitOfficial,fitTmp))
ssbplot(c(fitTmp,fitOfficial))

fbarplot(c(fitTmp,fitOfficial))

jj = jit(fitTmp,nojit = 20)
sim = simstudy(fit,nsim = 20)
plot(sim)
ret = retro(fit,year = 5)

plot(ret)

load("nssherring/fit_official.Rda")
fitOfficial = stockassessment:::refit(fit)
saveConf(fitOfficial$conf,file = "confOfficial.cfg")

stockassessment::predstdplot(fitTmp,fleet = 1, age =5:10)
corplot(fitTmp)













##Include error in data
aux = as.data.frame(dat$aux)

year = 2010
age = 5
fleet = 1

i = which(aux$year==year& aux$fleet == fleet& aux$age == age)
dat$logobs[i] = dat$logobs[i]+10

conf = loadConf(dat,"nssherring/suggestedConf/confTmp.cfg")
par = defpar(dat,conf)
fitWrong = sam.fit(dat,conf,par)
ssbplot(c(fitWrong,fitTmp))

conf = loadConf(dat,"nssherring/suggestedConf/confRobust.cfg")
par = defpar(dat,conf)
fitRobust = sam.fit(dat,conf,par)
ssbplot(c(fitFinal = fitTmp,fitRobust = fitRobust))



#fix model parameter
par$logSdLogN[2] = -2
map = list(logSdLogN = as.factor(c(0,NA)))
fitMap = sam.fit(dat,conf,par,map = map)
ssbplot(c(fitRobust,fitMap))
fitMap$pl

ii = which(dat$aux[,])
