# testing some survival analysis examples

library(survival)
dfData = read.csv(file='dataExternal/addicts.csv')
str(dfData)
## plot the kaplan-meier for relapse event
fit.km = survfit(Surv(survtime, status) ~ 1, conf.type='none', type='kaplan-meier', data=dfData)
plot(fit.km, mark=1)

## plot the kaplan-meier for relapse event
fit.km = survfit(Surv(survtime, status) ~ 1, type='kaplan-meier', data=dfData)
plot(fit.km, mark=1, conf.int=T)

## plot the kaplan-meier for relapse event stratified by covariates of choice
str(dfData)
fit.km = survfit(Surv(survtime, status) ~ clinic, type='kaplan-meier', data=dfData)
plot(fit.km, mark=1, conf.int=T, lty=1:2)
