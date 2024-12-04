library(stockassessment)
source("../utils/residualsTests.R")

cn<-read.ices("neahaddock/cn.dat")
cw<-read.ices("neahaddock/cw.dat")
dw<-read.ices("neahaddock/dw.dat")
lf<-read.ices("neahaddock/lf.dat")
lw<-read.ices("neahaddock/lw.dat")
mo<-read.ices("neahaddock/mo.dat")
nm<-read.ices("neahaddock/nm.dat")
pf<-read.ices("neahaddock/pf.dat")
pm<-read.ices("neahaddock/pm.dat")
sw<-read.ices("neahaddock/sw.dat")
surveys<-read.ices("neahaddock/survey.dat")

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

saveConf(confDef,file = "neahaddock/conf.cfg")

res = residuals(fit)
residDiagPlot(fit,resid = res)


######
conf01 = loadConf(dat,file = "neahaddock/suggestedConf/conf01Plus.cfg")
par01 = defpar(dat,conf01)
fit01 = sam.fit(dat,conf01,par01)

res01 = residuals(fit01)
residDiagPlot(fit01,resid = res01)

######
conf02 = loadConf(dat,file = "neahaddock/suggestedConf/conf02Qdens.cfg")
par02 = defpar(dat,conf02)
fit02 = sam.fit(dat,conf02,par02)
partable(fit02)
AIC(fit,fit01,fit02)
ssbplot(c(fit,fit01,fit02))

res02 = residuals(fit02)
residDiagPlot(fit02,resid = res02)

######
conf03 = loadConf(dat,file = "neahaddock/suggestedConf/conf03PredVar.cfg")
par03 = defpar(dat,conf03)
fit03 = sam.fit(dat,conf03,par03)
partable(fit03)
AIC(fit,fit01,fit02,fit03)
ssbplot(c(fit,fit01,fit02,fit03))

res03 = residuals(fit03)
residDiagPlot(fitTmp,resid = resTmp)

######
conf04 = loadConf(dat,file = "neahaddock/suggestedConf/conf04Cor.cfg")
par04 = defpar(dat,conf04)
fit04 = sam.fit(dat,conf04,par04)
partable(fit04)
AIC(fit,fit01,fit02,fit03,fit04)


res04 = residuals(fit04)
residDiagPlot(fit04,resid = res04)

jj = jit(fit04,nojit = 20)
sim = simstudy(fit04,nsim = 20)
ret = retro(fit04,year = 5)
proc = procres(fit04)

plot(ret)

dev.off()

######
confTmp = loadConf(dat,file = "neahaddock/suggestedConf/confTmp.cfg")
parTmp = defpar(dat,confTmp)
fitTmp = sam.fit(dat,confTmp,parTmp)
partable(fitTmp)
AIC(fit,fit01,fit02,fit03,fit04,fitTmp)

ssbplot(c(fit04,fitTmp),addCI = TRUE)
fbarplot(c(fit04,fitTmp),addCI = TRUE)
catchplot(fitTmp)
resTmp = residuals(fitTmp)
residDiagPlot(fitTmp,resid = resTmp)


fitWeb = stockassessment::fitfromweb("NEA_Had_2024")
fitWeb = stockassessment:::refit(fitWeb)

ssbplot(c(fitWeb = fitWeb, fit04 = fit04),addCI = TRUE)
fbarplot(c(fitWeb = fitWeb, fit04 = fit04),addCI = TRUE)
