library(stockassessment)
source("../utils/residualsTests.R")# Source script shared by Casper Berg for detecting significant structures in residuals

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

#saveConf(confDef,file = "nssherring/conf.cfg")
partable(fit)

#residDiagPlot(fit,resid = res)

conf01 = loadConf(dat,"nssherring/suggestedConf/conf01Fsta.cfg")
par01 = defpar(dat,conf01)
fit01 = sam.fit(dat,conf01,par01)
partable(fit01)

res01 = residuals(fit01)
residDiagPlot(fit01,resid = res01)


##########
conf02 = loadConf(dat,"nssherring/suggestedConf/conf02LogFpar.cfg")
par02 = defpar(dat,conf02)
fit02 = sam.fit(dat,conf02,par02)
partable(fit02)

AIC(fit,fit01,fit02)
res02 = residuals(fit02)
residDiagPlot(fit02,resid = res02)

##########
conf03 = loadConf(dat,"nssherring/suggestedConf/conf03PredVar.cfg")
par03 = defpar(dat,conf03)
fit03 = sam.fit(dat,conf03,par03)
partable(fit03)

AIC(fit,fit01,fit02,fit03)
res03 = residuals(fit03)
residDiagPlot(fit03,resid = res03)


##########
conf04 = loadConf(dat,"nssherring/suggestedConf/conf04ObsVar.cfg")
par04 = defpar(dat,conf04)
fit04 = sam.fit(dat,conf04,par04)
partable(fit04)

AIC(fit,fit01,fit02,fit03,fit04)
res04 = residuals(fit04)
residDiagPlot(fit04,resid = res04)


##########
conf05 = loadConf(dat,"nssherring/suggestedConf/conf05CorObs.cfg")
par05 = defpar(dat,conf05)
fit05 = sam.fit(dat,conf05,par05)
partable(fit05)

AIC(fit,fit01,fit02,fit03,fit04,fit05)
res05 = residuals(fit05)
residDiagPlot(fit05,resid = res05)


##########
confTmp = loadConf(dat,"nssherring/suggestedConf/confTmp.cfg")
parTmp = defpar(dat,confTmp)
fitTmp = sam.fit(dat,confTmp,parTmp)
partable(fitTmp)

#retTmp = retro(fitTmp,year = 5)

AIC(fit05,fitTmp)

AIC(fit,fit01,fit02,fit03,fit04,fit05,fitTmp)
resTmp = residuals(fitTmp)
residDiagPlot(fitTmp,resid = resTmp)




####
fitOfficial = readRDS("nssHerring/fitOfficial.Rda")
#saveConf(fitOfficial$conf,file = "nssHerring/confOfficial.cfg")

ssbplot(c(fitOfficial = fitOfficial,fitSuggested = fit05),addCI = TRUE)
fbarplot(c(fitOfficial = fitOfficial,fitSuggested = fit05),addCI = TRUE)
