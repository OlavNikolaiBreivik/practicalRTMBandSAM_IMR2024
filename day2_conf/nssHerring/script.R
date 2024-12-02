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

saveConf(confDef,file = "nssherring/conf.cfg")


#Usefull commands
AIC(fit) #Calcualte AIC
res = residuals(fit) #Caclulate osa-residuals
residDiagPlot(fit,resid = res) #Tests for patterns in residuals
jj = jit(fit) #Jitter analyuss
sim = simstudy(fit) #Simstudy
ret = retro(fit,year = 5) #Retrospective analysis
yearMat =  matrix(c(2023, 2024, 2024,2024, 
                    2022, 2023, 2023,2023, 
                    2021, 2022, 2022,2022,
                    2020, 2021, 2021,2021,
                    2019, 2020, 2020,2020), 
              nrow = 5, ncol = 4, byrow = TRUE) #year to remove in retor analysis
retManual = retro(fit,year = yearMat) #Retrospective analysis, manually provide years to remove from each fleet
