## 2024-05-10 discussion

library(survival)
library(timereg)


## pbc data from survival package
?pbc


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
summary(as.factor(dat$status))

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

summary(dat)

covs.hist = c("age", "logBilirubin", "logAlbumin", "logProtime")
par(mfrow=c(2,2))
for(i in 1:length(covs.hist)){
  hist(dat[,covs.hist[i]], xlab = covs.hist[i], main = "")
}







### Fit univairate Cox models --------------------------------------------------

# Fit a univariate model for sex
fit <-  coxph(Surv(time,status)~sex,
              data = dat)
fit
summary(fit)


# Fit a univariate model for age
fit<- coxph(Surv(time,status)~age,
            data = dat)
fit
summary(fit)




### Fit a multivariate Cox-PH model --------------------------------------------
fit_multi <- coxph(Surv(time,status)~sex+age+edema+logBilirubin+logAlbumin+logProtime,
             data = dat)
fit_multi
summary(fit_multi)


## Compute the R-square meansure
library(CoxR2)
coxr2fit = coxr2(fit_multi)
coxr2fit$rsq








## Check the PH assumption using cumulative martingale residuals
# Refit the Cox model with timereg:cox.aalen, which has cox regression as a special case.
library(timereg)
fit_multi2 <- cox.aalen(Surv(time,status)~prop(sex)+prop(age)+prop(edema)+
                          prop(logBilirubin)+prop(logAlbumin)+prop(logProtime),
                        data = dat)
fit_multi2

# p-values from test of PH assumption
summary(fit_multi2)

par(mfrow=c(2,3))
plot(fit_multi2, score = TRUE)








### Fit a stratified Cox-PH model with strata defined by edema
# Note that edema only takes 3 values: 0, 0.5, 1
summary(as.factor(dat$edema))

fit_strata <- coxph(Surv(time,status)~sex+age+strata(edema)+logBilirubin+logAlbumin+logProtime,
                    data = dat)
fit_strata    


## The cumulative baseline hazards for the stratified Cox model
bhaz <- basehaz(fit_strata, centered = FALSE)  # Get baseline hazards

# Plot
par(mfrow=c(1,1))
plot(bhaz$time[bhaz$strata=="edema=0"], bhaz$hazard[bhaz$strata=="edema=0"], 
     type = "s", col = 1, lty = 1, 
     xlab = "Time", ylab = "Cumulative Hazard", 
     ylim = range(bhaz$hazard),
     main = "Cumulative Baseline Hazard by Edema Strata")
lines(bhaz$time[bhaz$strata=="edema=0.5"], bhaz$hazard[bhaz$strata=="edema=0.5"], 
      type = "s", col = 2, lty = 2)
lines(bhaz$time[bhaz$strata=="edema=1"], bhaz$hazard[bhaz$strata=="edema=1"],
      type = "s", col = 3, lty = 3)
legend("topleft", legend = c("Edema = 0", "Edema = 0.5", "Edema = 1"),
       col = 1:3, lty = 1:3, bty = "n", title = "Strata")
# Here, type = "s" specifies that the line should be drawn as a step function, which is typical for survival plots.







### Fit a 3-piece piece-wise constant Cox model 

# pick cuts
qtl = quantile(dat$time[dat$status==1], c(1/3, 2/3))
qtl
qtl = as.numeric(qtl)
qtl

covs = c("sex", "age", "edema", "logBilirubin", "logAlbumin", "logProtime")
dat$sex = as.numeric(dat$sex== "f")   # 1: female, 0: male

# Create pseudo data with time-varying covariates Z*I_{(cut1, cut2]}(t)
dat2 = dat[dat$time> qtl[1], ]
dat3 = dat[dat$time> qtl[2], ]
covs0_1 = matrix(0, nrow = nrow(dat), ncol = length(covs))
covs0_2 = matrix(0, nrow = nrow(dat2), ncol = length(covs))
covs0_3 = matrix(0, nrow = nrow(dat3), ncol = length(covs))

covsname1 = paste(covs, "1", sep = "")
covsname2 = paste(covs, "2", sep = "")
covsname3 = paste(covs, "3", sep = "")

dat_pseudo1 = cbind(id = dat[,"id"], start = 0, stop = pmin(dat$time, qtl[1]), 
                   delta = as.numeric(dat$time<=qtl[1])*dat$status,
                   dat[, covs], covs0_1, covs0_1)
colnames(dat_pseudo1)[5:22] <- c(covsname1, covsname2, covsname3)

dat_pseudo2 = cbind(id = dat2[,"id"], start = qtl[1], stop = pmin(dat2$time, qtl[2]), 
                    delta = as.numeric(dat2$time<=qtl[2])*dat2$status,
                    covs0_2, dat2[, covs], covs0_2)
colnames(dat_pseudo2)[5:22] <- c(covsname1, covsname2, covsname3)

dat_pseudo3 = cbind(id = dat3[,"id"], start = qtl[2], stop = dat3$time, 
                    delta = dat3$status,
                    covs0_3, covs0_3, dat3[, covs])
colnames(dat_pseudo3)[5:22] <- c(covsname1, covsname2, covsname3)

dat_pseudo = rbind(dat_pseudo1, dat_pseudo2, dat_pseudo3)


# Change the order of roles - put the rows for same subjects together
order = order(as.numeric(dat_pseudo$id))
dat_psd = dat_pseudo[order, ]

# colnames(dat_psd)[5:22] <- c(covsname1, covsname2, covsname3)

dat_psd[1:10,]
head(dat)


# Fit the Cox model
fit_3piece = coxph(as.formula(paste("Surv(start, stop, delta)~",
                                    paste(c(covsname1, covsname2, covsname3), collapse = "+"))),
                   data = dat_psd, id = id)

fit_3piece
