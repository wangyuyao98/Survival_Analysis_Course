## 2024-05-13 discussion

## Vignette for Cox models with time dependent covariate and time dependent coeffcients
## https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf

library(survival)
library(timereg)

### Continue looking at the PBC data

data(pbc, package="survival")
summary(pbc)

# Consider the five explanatory variables that are found to be important by Fleming & Harrington (1991):
# age, edema, bilirubin, albumin, protime
# As a illustration, we also include an additional binary variable "sex"
covs = c("sex", "age", "edema", "bili", "albumin", "protime")
dat = pbc[, c("id", "time", "status", covs)]

# Remove subjects with missingness
dat = dat[complete.cases(dat), ]
dim(dat)

# View death as event and others as censoring
dat$status <- as.numeric(dat$status == 2)

# Change time to years
dat$time = dat$time/365.25


## Summary statistics
summary(dat)
# Histogram for the covariates
# covs.hist = c("age", "albumin", "protime")
covs.hist = c("age", "bili", "albumin", "protime")
par(mfrow=c(2,2))
for(i in 1:length(covs.hist)){
    hist(dat[,covs.hist[i]], xlab = covs.hist[i], main = "")
}

dat$logAlbumin = log(dat$albumin)
dat$logProtime = log(dat$protime)
dat$logBilirubin = log(dat$bili)


covs = c("sex", "age", "edema", "logBilirubin", "logAlbumin", "logProtime")


## Note that we can also use Surv() without changing the status variable to 0 and 1. 
## For example, we may do
fit = coxph(Surv(time, status==2)~age, data = pbc)
fit

# The above is the same as using the transformned status in dat
fit = coxph(Surv(time, status)~age, data = dat)
fit



#### In 2025-05-10 lab, we see one approach of fitting Cox model with piecewise constant effect by manually creating the pseudo data
# The following is an easier approach that creates pseudo data with survSplit() function
qtl = quantile(dat$time[dat$status==1], c(1/3, 2/3))
qtl = as.numeric(qtl)
qtl

dat_cut = survSplit(Surv(time,status)~., cut = qtl,
                     data = dat, episode = "tgroup")   
# This gives the same pseudo data as in the approach we saw last time.

# Fit the piecewise Cox model
# formula = as.formula(paste("Surv(tstart, time, status)~ strata(tgroup):", "(", paste(covs, collapse = "+"),")"))
# formula
fit_3p = coxph(Surv(tstart,time,status) ~ strata(tgroup):(sex+age+edema+logBilirubin+logAlbumin+logProtime),
               data = dat_cut, id = id)
fit_3p






### Test for proportional hazard assumption using Scheonfeld residuals
fit_cox = coxph(Surv(time,status)~sex+age+edema+logBilirubin+logAlbumin+logProtime,
                data = dat)
test.ph = cox.zph(fit_cox)
print(test.ph)
plot(test.ph)


# We can also do a graphical diagnostic using the function ggcoxzph() [in the survminer package], which produces, for each covariate, graphs of the scaled Schoenfeld residuals against the transformed time.
library(survminer)
ggcoxzph(test.ph)








