library(stockassessment)

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


#Usefull commands
AIC(fit)
res = residuals(fit)
residDiagPlot(fit,resid = res)
jj = jit(fit)
sim = simstudy(fit)
retro(fit(