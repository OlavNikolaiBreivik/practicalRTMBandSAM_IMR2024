fit = readRDS("forecast_ex/fitnssherring.Rda")

#Catch exactly 700 next year, and median fishing mortality is 0.2 for the two following years. Process for F is not propegated, asusmes same selectivity all years within each simulation.
fc = forecast(fit,catchval.exact = c(700,NA,NA),fval = c(NA,0.2,0.2), processNoiseF = FALSE,nosim = 1000)

fbarplot(fc)
ssbplot(fc)
fc

#Forecast 10 years forward in time given 700 catch this year and the same fishing mortality in all years.
fc = forecast(fit,catchval.exact = c(700,rep(NA,9)),fscale = c(NA,rep(1,9)), processNoiseF = FALSE)
fbarplot(fc)
ssbplot(fc)


