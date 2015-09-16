############################################################
### Code for module 5: the basic Nmix model for closed metapopulations
############################################################
# 5.3. Simulation of data and first analysis using unmarked

# Create a covariate called vegHt
nSites <- 100                # Also called M
set.seed(443)                # so that we all get the same values of vegHt
vegHt <- sort(runif(nSites, 1, 3)) # uniform from 1 to 3, sort for graph convenience

# Suppose that expected population size increases with vegHt
# The relationship is described by an intercept of -3 and
#    a slope parameter of 2 on the log scale
lambda <- exp(-3 + 2*vegHt)

# Now we go to 100 sites and observe the # of individuals (perfectly)
N <- rpois(nSites, lambda)

# We can fit a model that relates abundance to vegHt using the glm() function
#  with "family=Poisson":

summary(fm.glm1 <- glm(N ~ vegHt, family=poisson))

# Do some analysis of the results
plot(vegHt, N, xlab="Vegetation height", ylab="Abundance (N)")
glm1.est <- coef(fm.glm1)
plot(function(x) exp(-3 + 2*x), 1, 3, add=TRUE, lwd=3)
plot(function(x) exp(glm1.est[1] + glm1.est[2]*x), 1, 3, add=TRUE,
     lwd=3, col="blue")
legend(1, 20, c("Truth", "Estimate"), col=c("black", "blue"), lty=1,
       lwd=3)

nVisits <- 3
p <- 0.6
C <- matrix(NA, nSites, nVisits)
for(i in 1:nSites) {
    C[i,] <- rbinom(nVisits, N[i], p)
}

# Look at the data
cbind(N=N, C1=C[,1], C2=C[,2], C3=C[,3])

# Load library, format data and summarize
library(unmarked)
umf <- unmarkedFramePCount(y=C, siteCovs=as.data.frame(vegHt))
summary(umf)

# Fit a model and extract estimates
# Detection covariates follow first tilde, then come abundance covariates
# Also time model fit and monitor convergence
system.time(summary(fm.nmix1 <- pcount(~1 ~vegHt, data=umf, control=list(trace=TRUE, REPORT=1))))

# May consider other abundance models: NB = negative binomial, ZIP = zero-inflated Poisson (currently no others available)
# Compare with AIC
fm.nmix2<-  pcount(~1 ~vegHt, data=umf,mixture="NB", control=list(trace=TRUE, REPORT=1))
fm.nmix3<-  pcount(~1 ~vegHt, data=umf,mixture="ZIP", control=list(trace=TRUE, REPORT=1))
cbind(AIC.P=fm.nmix1@AIC, AIC.NB=fm.nmix2@AIC, AIC.ZIP=fm.nmix3@AIC)
        AIC.P   AIC.NB  AIC.ZIP
[1,] 870.3834 872.3845 872.3414


# Note, estimates of detection coefficients are on the logit-scale
# When there are no covariates, we can back-transform using:
beta1 <- coef(fm.nmix1) 
exp(beta1[3]) / (1+exp(beta1[3]))   # or
plogis(beta1[3])                    # or
backTransform(fm.nmix1, type="det") # estimate with SE

# When covariates are present we can do something like
plot(function(x) exp(beta1[1] + beta1[2]*x), 1, 3,
     xlab="vegetation height", ylab="Expected Abundance")

# Or suppose you want predictions for new values of vegHt, say 1.2 and 3.1
newdat <- data.frame(vegHt=c(1.2, 3.1))
predict(fm.nmix1, type="state", newdata=newdat)
> predict(fm.nmix1, type="state", newdata=newdat)
   Predicted         SE      lower      upper
1  0.6281663 0.09739854  0.4635472  0.8512464
2 24.5194200 2.70937154 19.7448241 30.4485852

ranef(fm.nmix1)

## 5.4. Analysis of Alder flycatcher point count data using unmarked
# PART 1: Set-up the data for analysis
# -------------------------- Format data ---------------------------------
# This a subset of point-count data from Chandler et al. (Auk 2009)
# alfl is Alder Flycatcher (Empidonax alnorum)

# Import data and check structure
#alfl.data <- read.csv("alfl05.csv", row.names=1)
alfl.data <- read.csv("http://sites.google.com/site/unmarkedinfo/home/webinars/2012-january/data/alfl05.csv?attredirects=0&d=1", row.names=1)
str(alfl.data)

# Pull out count matrix 
alfl.y <- alfl.data[,c("alfl1", "alfl2", "alfl3")]

# Standardize site-covariates
woody.mean <- mean(alfl.data$woody)
woody.sd <- sd(alfl.data$woody)
woody.z <- (alfl.data$woody-woody.mean)/woody.sd
struct.mean <- mean(alfl.data$struct)
struct.sd <- sd(alfl.data$struct)
struct.z <- (alfl.data$struct-struct.mean)/struct.sd

# Load and create unmarkedFrame 
# Note formatting of unmarked data frame
library(unmarked)
alfl.umf <- unmarkedFramePCount(y=alfl.y,
    siteCovs=data.frame(woody=woody.z, struct=struct.z),
    obsCovs=list(time=alfl.data[,c("time.1", "time.2", "time.3")],
                 date=alfl.data[,c("date.1", "date.2", "date.3")]))
summary(alfl.umf)

# Here's an easy way to standardize covariates after making the UMF
obsCovs(alfl.umf) <- scale(obsCovs(alfl.umf))
summary(alfl.umf)



# PART 2: Fit some models
# -------------------------- Model fitting  -----------------------------
(fm1 <-  pcount(~1 ~1, alfl.umf))
 backTransform(fm1, type="state")
 backTransform(fm1, type="det")
(fm2 <- pcount(~date+time ~1, alfl.umf))
(fm3 <- pcount(~date+time ~woody, alfl.umf))
(fm4 <- pcount(~date+time ~woody+struct, alfl.umf))
(fm5 <- pcount(~date+time ~1, alfl.umf,mixture="NB"))
(fm6 <- pcount(~date+time ~1, alfl.umf,mixture="ZIP"))
(fm7 <- pcount(~date+time ~woody,alfl.umf,mixture="ZIP"))
(fm8 <- pcount(~date+time ~struct,alfl.umf,mixture="ZIP"))
(fm9 <- pcount(~date+time ~woody+struct, alfl.umf,mixture="ZIP"))
(fm10<- pcount(~date+time ~woody+struct, alfl.umf,mixture="NB"))

# -------------------------- Model selection -----------------------------
# Put the fitted models in a "fitList"
fms <- fitList("lam(.)p(.)"                    = fm1,
               "lam(.)p(date+time)"            = fm2,
               "lam(woody)p(date+time)"        = fm3,
               "lam(woody+struct)p(date+time)" = fm4,
               "lam(.)p(date+time)NB"          = fm5,
               "lam(.)p(date+time)ZIP"         = fm6,
               "lam(woody)p(date+time)ZIP"     = fm7,
               "lam(struct)p(date+time)ZIP"    = fm8,
               "lam(woody+struct)p(date+time)ZIP"=fm9,
               "lam(woody+struct)p(date+time)NB" =fm10)

# Rank them by AIC
(ms <- modSel(fms))

# Table with everything you could possibly need
coef(ms)
toExport <- as(ms, "data.frame")
str(toExport)               # See what?s in there



# PART 3: Do some analysis of the results
# ---------------------------- Prediction --------------------------------
# Expected detection probability as function of time of day
# We standardized "time", so we predict over range of values on that scale
# We must fix "date" at some arbitrary value (let's use the mean)
newData1 <- data.frame(time=seq(-2.08, 1.86, by=0.1), date=0)
E.p <- predict(fm4, type="det", newdata=newData1, appendData=TRUE)
head(E.p)

# Plot it
# Expected detection over time of day
par(mfrow=c(1,2))
plot(Predicted ~ time, E.p, type="l", ylim=c(0,1), xlab="time of day (standardized)", ylab="Expected detection probability", col = "blue", lwd = 3)
lines(lower ~ time, E.p, type="l", col=gray(0.5))
lines(upper ~ time, E.p, type="l", col=gray(0.5))

# Expected abundance over range of "woody"
newData2 <- data.frame(woody=seq(-1.6, 2.38,,50),struct=seq(-1.8,3.2,,50))
E.N <- predict(fm4, type="state", newdata=newData2, appendData=TRUE)
head(E.N)

# Plot predictions with 95% CI
plot(Predicted ~ woody, E.N, type="l", ylim=c(-.1,max(E.N$Predicted)), xlab="woody vegetation (standardized)", ylab="Expected abundance, E[N]", col = "blue", lwd = 3)
lines(lower ~ woody, E.N, type="l", col=gray(0.5))
lines(upper ~ woody, E.N, type="l", col=gray(0.5))

# Plot it again, but this time convert the x-axis back to original scale
par(mfrow=c(1,1))
plot(Predicted ~ woody, E.N, type="l", ylim=c(-.1,max(E.N$Predicted)),
     xlab="Percent cover - woody vegetation",
     ylab="Expected abundance, E[N]",
     xaxt="n")
xticks <- -1:2
xlabs <- xticks*woody.sd + woody.mean
axis(1, at=xticks, labels=round(xlabs, 1))
lines(lower ~ woody, E.N, type="l", col=gray(0.5))
lines(upper ~ woody, E.N, type="l", col=gray(0.5))

# Get predictions for covariates in the observation model on natural scale
# More step-by-step approach
# Look at range, mean and sd of covariate
(range.time <- range(alfl.data[,7:9]))
(mean.time <- mean(as.matrix(alfl.data[,7:9])))
(sd.time <- sd(c(as.matrix(alfl.data[,7:9]))))
(range.date <- range(alfl.data[,10:12]))
(mean.date <- mean(as.matrix(alfl.data[,10:12])))
(sd.date <- sd(c(as.matrix(alfl.data[,10:12]))))

# Create new covariate for prediction and scale identically
original.time.pred <- seq(5,10,,100)
original.date.pred <- seq(1,54,,100)
time.pred <- (original.time.pred - mean.time) / sd.time
date.pred <- (original.date.pred - mean.date) / sd.date

# Compute predictions for both covariates, keeping other constant
newData1 <- data.frame(time=time.pred, date=0)
Ep1 <- predict(fm4, type="det", newdata=newData1, appendData=TRUE)
head(Ep1)
newData2 <- data.frame(time=0, date=date.pred)
Ep2 <- predict(fm4, type="det", newdata=newData2, appendData=TRUE)
head(Ep2)

# Plot against covariate on natural scale
par(mfrow=c(1,2))
plot(Ep1$Predicted ~ original.time.pred, type="l", ylim=c(0,1), xlab="Time of day", ylab="Expected detection probability", main = "Effect of time of day", col = "blue", lwd = 3)
lines(Ep1$lower ~ original.time.pred, type="l", col=gray(0.5))
lines(Ep1$upper ~ original.time.pred, type="l", col=gray(0.5))

plot(Ep2$Predicted ~ original.date.pred, type="l", ylim=c(0,1), xlab="Date", ylab="Expected detection probability", main = "Effect of date", col = "blue", lwd = 3)
lines(Ep2$lower ~ original.date.pred, type="l", col=gray(0.5))
lines(Ep2$upper ~ original.date.pred, type="l", col=gray(0.5))

# ---------------------------- Goddess of fit -------------------------------
### Here's an example of a bootstrap GoF analysis.
### Best model is in "fm4" object

# Function returning three fit-statistics.
fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2)                                 # Sums of squares
    chisq <- sum((observed - expected)^2 / expected)     # Chisquared
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2) # Freeman-Tukey
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
    }

(pb <- parboot(fm4, fitstats, nsim=100, report=1))
print(pb)
plot(pb)

## Now let?s bootstrap a summary statistic 
## This is not too meaningful right now but we will do a similar thing 
## later in a more relevant context

# Total population size (derived parameter) over 50 surveyed sites
Nhat <- function(fm) {
    N <- sum(predict(fm, type="state")$Predicted, na.rm=TRUE)
    }
    
(pb.N <- parboot(fm4, Nhat, nsim=25, report=5))
plot(pb.N)

# Here's an example of model-averaging predictions
# See pg 150, section 4.2.1, of Burnham and Anderson (2002)
# This might be worthwhile since fm3 and fm4 had similar support
newData3 <- data.frame(woody=seq(-1.6, 2.38,,50), struct=seq(-1.8,3.2,,50))

## averages over _all_ models in the fit list "fms":
E.N.bar <- predict(fms, type="state", newdata=newData3, appendData=TRUE)
head(E.N.bar)

# Plot it
plot(Predicted ~ woody, E.N.bar, type="l", ylim=c(-0.1, max(E.N$Predicted)),
     xlab="Percent cover - woody vegetation",
     ylab="Expected abundance, E[N]",
     xaxt="n")
xticks <- -1:2
xlabs <- xticks*woody.sd + woody.mean
axis(1, at=xticks, labels=round(xlabs, 1))
lines(lower ~ woody, E.N.bar, type="l", col=gray(0.5))
lines(upper ~ woody, E.N.bar, type="l", col=gray(0.5))



## 5.5. Analysis of Swiss willow tits with unmarked

# PART 1: Set-up of analysis
# Read in some data from the Swiss MHB survey (Swiss equivalent of BBS)
mhbdata <- read.csv("http://sites.google.com/site/unmarkedinfo/home/webinars/2012-january/data/wtmatrix.csv?attredirects=0&d=1")
#mhbdata<-read.csv("wtmatrix.csv")
mhbdata[1:10,]

library("unmarked")

mhb.y<-mhbdata[,c("c.1","c.2","c.3")]
mhbdata[,"length"]<-1/mhbdata[,"length"]
mhb.umf<-unmarkedFramePCount(y=mhb.y,
siteCovs=data.frame(elev=mhbdata[,"elev"],forest=mhbdata[,"forest"],length=mhbdata[,"length"]),
obsCovs=list(duration=mhbdata[,c("dur.1","dur.2","dur.3")],
             day = mhbdata[,c("day.1","day.2","day.3")])  )
# this is extremely handy:
obsCovs(mhb.umf)<- scale(obsCovs(mhb.umf))
#
# NOTE: Do not standardize 1/length because we are using 1/length
# for a specific reason
siteCovs(mhb.umf)$forest<-scale(siteCovs(mhb.umf)$forest)
siteCovs(mhb.umf)$elev<-scale(siteCovs(mhb.umf)$elev)
str(mhb.umf)


## PART 2: Fit some models
# ------ Fit Poisson N-mixture model to MHB data and do model selection -----#
## Do simplistic two-step model selection with p first
fm01 <- pcount(~day ~1, mhb.umf)
fm02 <- pcount(~day+I(day^2) ~1, mhb.umf)
fm1<- pcount(~1 ~1, mhb.umf)

# Put three fitted models in a "fitList" and rank them by AIC
fms <- fitList("lam(.)p(.)"             = fm1,
               "lam(.)p(day)"           = fm01,
               "lam(.)p(day+day2)"      = fm02)
(ms <- modSel(fms))
# So no evidence for seasonal effect on detection within simplest model for lambda

# Go on modeling abundance part of model with p=constant
fm2 <- pcount(~1 ~elev,mhb.umf, control=list(trace=TRUE, REPORT=1))
fm3<- pcount(~1 ~forest,mhb.umf, control=list(trace=TRUE, REPORT=1))
fm4<- pcount(~1 ~length,mhb.umf, control=list(trace=TRUE, REPORT=1))
fm5<- pcount(~1 ~forest+elev,mhb.umf, control=list(trace=TRUE, REPORT=1))
fm6<- pcount(~1 ~forest+length,mhb.umf, control=list(trace=TRUE, REPORT=1))
fm7<- pcount(~1 ~elev+length,mhb.umf, control=list(trace=TRUE, REPORT=1))
fm8<- pcount(~1 ~forest+elev+length,mhb.umf, control=list(trace=TRUE, REPORT=1))
fm9<- pcount(~1 ~elev + I(elev^2),mhb.umf, control=list(trace=TRUE, REPORT=1))
fm10<- pcount(~1 ~forest+elev+I(elev^2) + length,mhb.umf, control=list(trace=TRUE, REPORT=1))
fm11<- pcount(~1 ~forest+elev+I(elev^2)+I(elev^3) + length,mhb.umf, control=list(trace=TRUE, REPORT=1))

# Put fitted models in a "fitList" and rank them by AIC
mspart1<- fitList(
"lam(.)p(.)"                         = fm1,
"lam(elev)p(.)"                      = fm2,
"lam(forest)p(.)"                    = fm3,
"lam(length)p(.)"                    = fm4,
"lam(forest+elev)p(.)"               = fm5,
"lam(forest+length)p(.)"             = fm6,
"lam(elev+length)p(.)"               = fm7,
"lam(forest+elev+length)p(.)"        = fm8,
"lam(elev + elev^2)p(.)"             = fm9,
"lam(forest+elev+elev^2+length)p(.)" = fm10,
"lam(forest+elev+elev^2+elev^3+length)p(.)" = fm11)
(ms1 <- modSel(mspart1))

print(coef(ms1), dig = 2)

# Revisit date effects within current best model and compare with AIC
system.time(fm12 <- pcount(~day ~forest+elev+I(elev^2) +I(elev^3) + length,
     mhb.umf,control=list(trace=TRUE, REPORT=1)))
system.time(fm13 <- pcount(~day + I(day^2)~forest+elev+I(elev^2)+ I(elev^3) +
     length, mhb.umf, control=list(trace=TRUE, REPORT=1)))
mspart2<- fitList(
"lam(forest+elev+elev^2+elev^3+length)p(.)"           = fm10,
"lam(forest+elev+elev^2+elev^3+length)p(day)"         = fm12,
"lam(forest+elev+elev^2+elev^3+length)p(day+day2)"    = fm13)
(ms2 <- modSel(mspart2))
# Now good evidence for seasonal effects on p

# Table with everything you could possibly need
coef(ms2)
(toExport <- as(ms2, "data.frame"))



### PART 3: Do some analysis of the results
# --------- Goodness of fit of Poisson models -----------#
## Does this model fit worth a darn?
# Function returning three fit-statistics.
# NOTE: na.rm=TRUE !!!!!!
fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2,na.rm=TRUE)
    chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
    }

# (pb.mhb <- parboot(fm13, fitstats, nsim=100, report=1)) # this takes a while
(pb.mhb <- parboot(fm13, fitstats, nsim=10, report=1)) # cheaper version
plot(pb.mhb)
# Ouch ! Model does not fit!


# What should we do ? 
# Three choices:
# 1. We can expand this model or tinker with components of it (NB, ZIP)
#   [takes a long time to run these models]
# 2. We can seek out an implausible data-generating model that fits.
#   [might satisfy referee to have good p-value]
# 3. We can proceed anyway.


# Try choice 1:
# --- Fit NegBin N-mixture models to MHB data and do model selection ----#
# Try same models with NegBin N: accounts for extra-Poisson dispersion
system.time(fm20 <- pcount(~1 ~elev, mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))
system.time(fm21<- pcount(~1 ~forest, mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))
system.time(fm22<- pcount(~1 ~length, mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))
system.time(fm23<- pcount(~1 ~forest+elev, mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))
system.time(fm24<- pcount(~1 ~forest+length, mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))
system.time(fm25<- pcount(~1 ~elev+length, mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))
system.time(fm26<- pcount(~1 ~forest+elev+length, mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))
system.time(fm27<- pcount(~1 ~elev + I(elev^2), mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))
system.time(fm28<- pcount(~1 ~forest+elev+I(elev^2) + length, mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))
system.time(fm29<- pcount(~1 ~forest+elev+I(elev^2)+I(elev^3) + length, mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))
system.time(fm30<- pcount(~day ~forest+elev+I(elev^2)+I(elev^3) + length, mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))
system.time(fm31<- pcount(~day + I(day^2) ~forest+elev+I(elev^2)+I(elev^3) + length, mixture = "NB", mhb.umf, control=list(trace=TRUE, REPORT=1)))


# Put fitted models in a "fitList" and rank them by AIC
# Add best Poisson model for comparison
mspart3<- fitList(
"lam(P, forest+elev+elev^2+elev^3+length)p(day+day^2)" = fm13,
"lam(NB, elev)p(.)"                      = fm20,
"lam(NB, forest)p(.)"                    = fm21,
"lam(NB, length)p(.)"                    = fm22,
"lam(NB, forest+elev)p(.)"               = fm23,
"lam(NB, forest+length)p(.)"             = fm24,
"lam(NB, elev+length)p(.)"               = fm25,
"lam(NB, forest+elev+length)p(.)"        = fm26,
"lam(NB, elev + elev^2)p(.)"             = fm27,
"lam(NB, forest+elev+elev^2+length)p(.)" = fm28,
"lam(NB, forest+elev+elev^2+elev^3+length)p(.)" = fm29,
"lam(NB, forest+elev+elev^2+elev^3+length)p(day)" = fm30,
"lam(NB, forest+elev+elev^2+elev^3+length)p(day+day^2)" = fm31)
(ms3 <- modSel(mspart3))

# --- Do parboot GOF test for current AIC best NegBin model -----------#
 (pb.mhb <- parboot(fm31, fitstats, nsim=10, report=1)) # cheap version
plot(pb.mhb)
# Looks much better now. Use this model for inference about willow tit distribution and abundance in Switzerland


# ---------------------- Do 1D predictions ----------------------------#
## Now we want to use these parameter estimates to do a couple things
## Predict response of E[N] vs. elevation -- find optimal elevation
## make an abundance map over Switzerland 

range(mhbdata[,"elev"])
> range(mhbdata[,"elev"])
[1]  250 2750
(elev.mean<- attr(siteCovs(mhb.umf)$elev,"scaled:center"))
(elev.sd <- attr(siteCovs(mhb.umf)$elev,"scaled:scale"))

# remember length = 0 is saturation sampling because length = 1/L
# Create covariate values for prediction and then scale them identically
original.pred.elev <- seq(250,2750,,1000)
pred.elev <- (original.pred.elev - elev.mean) / elev.sd
newData<- data.frame(elev=pred.elev, forest=0, length=0)
pred<-predict(fm31, type="state", newdata=newData, appendData=TRUE)
head(pred)

plot(Predicted ~ original.pred.elev, pred,type="l",xlab="Elevation (not standardized)", ylab="Expected # territories",ylim=c(0,100), lwd = 2)
lines(lower ~ original.pred.elev, pred,type="l",col="red", lwd = 2)
lines(upper ~ original.pred.elev, pred,type="l",col="red", lwd = 2)
# Note considerable uncertainty with NB state model

# what is primo elevation for the willow tit?

# For quadratic response could use calculus
#   quadratic response:
#      y = a + b*x + c*x2
#   differentiate and set to 0:
#      dy/dx = b + 2*c*x = 0
#   solve
#      xopt = -b/(2*c)

# Empirical answer for any model of covariate
(original.pred.elev[pred$Predicted == max(pred$Predicted)])
[1] 1891.642


### PART 3c:
### SPATIAL ANALYSIS/PREDICTION
### Now lets make some really cool spatial predictions (aka MAPS)

landscape <- read.csv("http://sites.google.com/site/unmarkedinfo/home/webinars/2012-january/data/Swiss_landscape.csv?attredirects=0&d=1")
# landscape<-read.csv("Swiss_landscape.csv")
head(landscape)
## Note integer coordinates - row/column ids

gelev<- landscape[,"medel"]   # median elevation of quadrat
gforest<-landscape[,"forest"]
grid<-landscape[,c("x","y")]

# lets plot these variables to see how they look
#
# two options: (1) use my simpleton spatial.plot function
#              (2) stuff the data into a matrix and use image()
#
# grab utility functions including spatial.plot

source("http://sites.google.com/site/unmarkedinfo/home/webinars/2012-january/data/utils.R?attredirects=0&d=1")

# par(mar=c(3,3,3,5),mfrow=c(2,1))
par(mar=c(3,3,3,6))
spatial.plot(grid,gelev)
text(500, 290, labels = "Elevation")
spatial.plot(grid,gforest)
text(500, 290, labels = "Forest cover")
# this is cool


# ------------------ Do 2D predictions (maps) ----------------------------#
library(sp)
# Get coordinates and rescale from km to m
coordCH <- matrix(cbind(landscape$x + 0.5, landscape$y + 0.5), ncol = 2)
xcor <- coordCH[,1] * 1000
ycor <- coordCH[,2] * 1000

# Get predictions for each pixel of Switzerland
# First scale values of elevation analogously to analysis
pelev <- (landscape$medel-elev.mean)/elev.sd
forest.mean <- mean(mhbdata[,"forest"])
forest.sd <- sd(mhbdata[,"forest"])
pforest <- (landscape$forest-forest.mean)/forest.sd
### new<- data.frame(elev=pelev,forest=pforest,length=0)
## pred<-predict(fm31,type="state",newdata=new,appendData=TRUE)
###
### this would destroy your computer -> bug in unmarked
###
### instead we have to do this the old-fashioned way:
### look at col names to figure order of parameters

par(mfrow=c(1,1))
betavec<-coef(fm31)[1:5]
Xg<-cbind(rep(1,length(pelev)),pforest,pelev,pelev^2, pelev^3)
pred<-exp(Xg%*%(betavec))

# Define a new dataframe with coordinates and outcome to be plotted 
PARAM <- data.frame(x = xcor, y = ycor, z = pred)

# Convert the dataframe first into a SpatialPointsDataFrame and then into a SpatialPixelsDataFrame
coordinates(PARAM)<- ~x+y
gridded(PARAM) <- TRUE

# Plot the map using custom color palette
# mapPalette <- colorRampPalette(c("grey", "lightgreen", "darkgreen"))
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
spplot(PARAM, col.regions = mapPalette(100), main = "Expected willow tit density")
# Perhaps too high ? check with literature


# Next get total population size of breeding territories (derived parameter). 
N<- sum(pred)
print(N)

# Bootstrap the SE
# note: not using predict() here but doing calculation by hand
Nhat <- function(fm) {
   betavec<-coef(fm)[1:5]
   Xg<-cbind(rep(1,length(pelev)),pforest,pelev,pelev^2, pelev^3)
   pred<-exp(Xg%*%(betavec))
   N<-sum(pred)
   N
   }

estimate.of.territories<-Nhat(fm31)
(pb.N <- parboot(fm31, Nhat, nsim=100, report=1))  # Would have to do many more times
plot(pb.N)
bs.sample <- pb.N@t.star
summary(bs.sample)
quantile(bs.sample, prob = c(0.025, 0.975))
> quantile(bs.sample, prob = c(0.025, 0.975))
     2.5%     97.5% 
 132110.9 1193516.1                # 95% confidence interval
# Much uncertainty in estimate !

     
############################################################################################
## 5.6. Analysis of the N-mixture model with BUGS/JAGS and introduction to Bayesian p-values
# Bundle data
win.data <- list(C = C, R = nrow(C), T = ncol(C), vegHt = vegHt)


# Specify model in BUGS language
sink("Nmix.txt")
cat("
model {

# Priors
alpha0 ~ dunif(-10, 10)
alpha1 ~ dunif(-10, 10)
p ~ dunif(0, 1)

# Likelihood
# Ecological model for true abundance
for (i in 1:R){
   N[i] ~ dpois(lambda[i])
   log(lambda[i]) <- alpha0 + alpha1 * vegHt[i]

   # Observation model for replicated counts
   for (j in 1:T){
      C[i,j] ~ dbin(p, N[i])

         # Assess model fit using Chi-squared discrepancy
         # Compute fit statistic E for observed data
         eval[i,j] <- p * N[i]   	# Expected values
         E[i,j] <- pow((C[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
         # Generate replicate data and compute fit stats for them
         C.new[i,j] ~ dbin(p, N[i])
         E.new[i,j] <- pow((C.new[i,j] - eval[i,j]),2) / (eval[i,j] + 0.5)
   }
}

# Derived quantities
lp <- logit(p)           # logit-scale detection
totalN <- sum(N[])       # Total abundance over all R sites
fit <- sum(E[,])         # Fit stats actual data set
fit.new <- sum(E.new[,]) # Fit stats ?ideal? data set
}
",fill = TRUE)
sink()


# Define function to generate random initial values
Nst <- apply(C, 1, max) + 1	# Can be important to give good inits for N
inits <- function() list(N = Nst, alpha0 = runif(1, -1, 1), alpha1 = runif(1, -1, 1), p = runif(1))

# Parameters monitored
params <- c("alpha0", "alpha1", "p", "lp", "totalN", "fit", "fit.new")

# MCMC settings
ni <- 10000
nt <- 2
nb <- 1000
nc <- 3

# Call WinBUGS from R (ART 0.7 min)
library(R2WinBUGS)
out1 <- bugs(win.data, inits, params, "Nmix.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, working.directory = getwd())

# Summarize posteriors
print(out1, 3)

# Compare with MLEs from unmarked
summary(fm.nmix1 <- pcount(~1 ~vegHt, data=umf))

# Call JAGS from R (ART 1.26 min)
library("R2jags")		# requires rjags
system.time(out2 <- jags(win.data, inits, params, "Nmix.txt", n.chains = nc,
   n.thin = nt, n.iter = ni, n.burnin = nb) )
traceplot(out2)

# Summarize posteriors
print(out2, dig = 3)

par(mar = c(5,5,4,5), cex.lab = 1.5)
plot(vegHt, N, xlab="Vegetation height", ylab="Abundance (N)", las = 1)
glm1.est <- coef(fm.glm1)
plot(function(x) exp(-3 + 2*x), 1, 3, add=TRUE, lwd=3)
lines(vegHt, predict(fm.nmix1, "state")[,1], lwd=3, col="blue")
lines(vegHt, exp(out1$mean$alpha0 + out1$mean$alpha1 * vegHt), lwd=3, col="green", lty = "dashed")
lines(vegHt, exp(out2$BUGSoutput$mean$alpha0 + out2$BUGSoutput$mean$alpha1 * vegHt), lwd=3, col="red", lty = "dotted")
legend(1, max(N), c("Truth", "unmarked", "WinBUGS", "JAGS"), col=c("black", "blue", "green", "red"), lty=c(1, 1, 2, 3), lwd=3)

plot(out1$sims.list$fit, out1$sims.list$fit.new, xlab="Actual data set", ylab="?Perfect? data sets", las = 1)
abline(0,1, col = "red", lwd = 3)

(bpv <- mean(out1$sims.list$fit.new > out1$sims.list$fit))
# (bpv <- mean(out1$sims.list$fit.new > out1$sims.list$fit))
[1] 0.2846667


######################################################################################
tmp[i] <- step(10-N[i])
sum.critical <- sum(tmp[])	# Number of pops with critical size

plot(table(out1$sims.list$sum.critical), xlab="Number of threshold populations", ylab="Frequency")
abline(v = 87.9, col = "red", lwd = 5)

(metapop.extinction.risk <- mean(out1$sims.list$sum.critical>87))

## 5.10. Exercises

# Solution A:
range(mhbdata[,12:14], na.rm = TRUE)
day.mean <- mean(as.matrix(mhbdata[,12:14]), na.rm = TRUE)
day.sd <- sd(c(as.matrix(mhbdata[,12:14])), na.rm = TRUE)
original.pred.day <- 15:110
pred.day <- (original.pred.day - day.mean) / day.sd
new<- data.frame(day=pred.day)
pred<-predict(fm31,type="det",newdata=new,appendData=TRUE)
head(pred)

plot(Predicted ~ original.pred.day, pred,type="l",xlab="Date (1 = 1 April)", ylab="Expected detection prob",ylim=c(0,1), lwd = 2)
lines(lower ~ original.pred.day, pred,type="l",col="red", lwd = 2)
lines(upper ~ original.pred.day, pred,type="l",col="red", lwd = 2)