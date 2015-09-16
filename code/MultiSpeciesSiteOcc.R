


#### code from: http://www.esajournals.org/doi/abs/10.1890/0012-9658%282006%2987%5B842%3AESRAAB%5D2.0.CO%3B2
#### full article here: http://www.uvm.edu/rsenr/vtcfwru/spreadsheets/occupancy/Occupancy%20Exercises/Exercise15/Dorazio_et_al_2006.pdf
#### example here http://www.esapubs.org/Archive/ecol/E087/050/suppl-1.htm 


# library(R2WinBUGS)

# The model
sink("MultiSpeciesSiteOccModel")
cat("
model {

omega ~ dunif(0,1)
    
    psi.mean ~ dunif(0,1)
    alpha <- log(psi.mean) - log(1-psi.mean)
    
    theta.mean ~ dunif(0,1)
    beta <- log(theta.mean) - log(1-theta.mean)
    
    tau.u ~ dgamma(0.1,0.1)
    tau.v ~ dgamma(0.1,0.1)
    rho ~ dunif(-1,1)
    var.eta <- tau.v/(1.-pow(rho,2))
    
    sigma.u <- 1/sqrt(tau.u)
    sigma.v <- 1/sqrt(tau.v)
    
    
    for (i in 1:(n+nzeroes)) {
    w[i] ~ dbin(omega, 1)
    phi[i] ~ dnorm(alpha, tau.u)
    
    mu.eta[i] <- beta + (rho*sigma.v/sigma.u)*(phi[i] - alpha)
    eta[i] ~ dnorm(mu.eta[i], var.eta)
    
    
    logit(psi[i]) <- phi[i]
    logit(theta[i]) <- eta[i]
    
    mu.psi[i] <- psi[i]*w[i]
    for (j in 1:J) {
    Z[i,j] ~ dbin(mu.psi[i], 1)
    mu.theta[i,j] <- theta[i]*Z[i,j]
    X[i,j] ~ dbin(mu.theta[i,j], K)
    }
    }
    
    n0 <- sum(w[(n+1):(n+nzeroes)])
    N <- n + n0
    }
    ", fill=TRUE)
sink()

############################

MultiSpeciesSiteOcc <- function(nrepls, X) {
  
  start.time = Sys.time()
  
  # augment data matrix with an arbitrarily large number of zero row vectors
  nzeroes = 100
  n = dim(X)[1]
  nsites = dim(X)[2]
  Xaug = rbind(X, matrix(0, nrow=nzeroes, ncol=nsites))
  
  # create arguments for bugs()
  sp.data = list(n=n, nzeroes=nzeroes, J=nsites, K=nrepls, X=Xaug)
  
  sp.params = list('alpha', 'beta', 'rho', 'sigma.u', 'sigma.v', 'omega', 'N')
  
  sp.inits = function() {
    omegaGuess = runif(1, n/(n+nzeroes), 1)
    psi.meanGuess = runif(1, .25,1)
    theta.meanGuess = runif(1, .25,1)
    rhoGuess = runif(1, 0,1)
    sigma.uGuess = 1
    sigma.vGuess = 1
    list(omega=omegaGuess, psi.mean=psi.meanGuess, theta.mean=theta.meanGuess, tau.u=1/(sigma.uGuess^2), tau.v=1/(sigma.vGuess^2), rho=rhoGuess,
         w=c(rep(1, n), rbinom(nzeroes, size=1, prob=omegaGuess)),
         phi=rnorm(n+nzeroes, log(psi.meanGuess/(1.-psi.meanGuess)), sigma.uGuess),
         eta=rnorm(n+nzeroes, log(theta.meanGuess/(1.-theta.meanGuess)), sigma.vGuess),
         Z = matrix(rbinom((n+nzeroes)*nsites, size=1, prob=psi.meanGuess), nrow=(n+nzeroes))
    )
  }
  
  # fit model to data using WinBUGS code calling openbugs
  library(R2WinBUGS)
    
  fit = bugs(sp.data, sp.inits, sp.params,
             model.file='MultiSpeciesSiteOccModel',# debug=F, 
             n.chains=5, n.iter=55000, n.burnin=500, n.thin=50,
             program="openbugs")
  

# # OpenBugs
#   library ("BRugs")
#   fit = BRugsFit(data=sp.data, inits=sp.inits, parametersToSave=sp.params,
#              modelFile='MultiSpeciesSiteOccModel', coda = TRUE,
#              DIC = TRUE, 
#              numChains=4, nIter=5500, nBurnin=500, nThin=50)
#              

# # Parameters monitored
# sp.params = c('alpha', 'beta', 'rho', 'sigma.u', 'sigma.v', 'omega', 'N')

# 
# 
# # MCMC settings
# ni <- 5500
# nt <- 50
# nb <- 500
# nc <- 3
# 
#    # library("rjags")
#    library(R2jags)
#    fit <- jags(data=sp.data,  parameters.to.save = sp.params, 
#                model.file = "MultiSpeciesSiteOccModel", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
#    

  
  end.time = Sys.time()
  elapsed.time = difftime(end.time, start.time, units='mins')
  cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes\n', sep=''))
  
  list(fit=fit, data=sp.data, X=X)
}

##############################
# Cum Num Species Present ###
###############################

CumNumSpeciesPresent = function(nsites, alpha, sigmaU, N) {
  
  # Computes a sample of the posterior-predictive distribution of the (cumulative) number of species present at nsites.
  
  # compute posterior predictions of species occurrence probabilities
  ndraws = length(alpha)
  Nmax = max(N)
  logitPsi = matrix(NA, nrow=ndraws, ncol=Nmax)
  psi = logitPsi
  for (i in 1:ndraws) {
    logitPsi[i,1:N[i]] = rnorm(N[i], mean=alpha[i], sd=sigmaU[i])
    psi[i, 1:N[i]] = 1/(1+exp(-logitPsi[i,1:N[i]]))
  }
  
  # compute posterior predictions of species presence at each site
  z = array(NA, dim=c(ndraws, Nmax, nsites))
  for (i in 1:ndraws) {
    for (j in 1:N[i]) {
      z[i,j, ] = rbinom(nsites, size=1, prob=psi[i,j])
    }
  }
  
  # compute posterior predictions of cumulative number of species present
  M = matrix(NA, nrow=ndraws, ncol=nsites)
  for (i in 1:ndraws) {
    for (j in 1:nsites) {
      zsum = rep(NA, N[i])
      if (j>1) {
        zsum = apply(z[i, 1:N[i], 1:j], 1, sum)
      }
      else {
        zsum = z[i, 1:N[i], 1]
      }
      M[i,j] = sum(zsum>0)
    }
  }
  
  # compute summary stats for plotting
  nSpeciesPresent = matrix(NA, nrow=3, ncol=nsites)
  for (j in 1:nsites) {
    x = M[,j]
    nSpeciesPresent[1, j] = mean(x)
    nSpeciesPresent[2:3, j] = quantile(x, probs=c(.05, .95))
  }
  
  # plot results
  ylimits = c(min(nSpeciesPresent[2,]), max(nSpeciesPresent[3,]))
  plot(1:nsites, nSpeciesPresent[1,], pch=16, ylim=ylimits, type='b',
       xlab='Number of sampled plots', ylab='Number of species', las=1, cex.axis=1.2, cex.lab=1.5, cex=1.5)
  segments(1:nsites, nSpeciesPresent[2,], 1:nsites, nSpeciesPresent[3,])
  
  list(meanAndquantiles=nSpeciesPresent, summaryStats=summary(M))
}




# 
# ##########################
# #  apply model ###########
# #########################
# 
# X1 = as.matrix(read.csv("data/PS1.csv"))
# nrepls = 12
# algas = MultiSpeciesSiteOcc(nrepls, X1)
# 
# algas$fit
# library(mcmcplots)
# mcmcplot(algas$fit)
# 
# alpha.post = algas$fit$sims.matrix[,"alpha"]
# sigmaU.post = algas$fit$sims.matrix[,"sigma.u"]
# N.post = algas$fit$sims.matrix[,"N"]
# 
# nsites = 35
# CumNumSpeciesPresent(nsites, alpha.post, sigmaU.post, N.post)
# 
# 
# 
# #########################################
# ######## from Bayesian Pop Anal. Kery book
# ########################################
# 
# # 6.3. Analysis of a real data set: model Mtbh for species richness estimation
# # Read in data and look at them
# # p610 <- read.table("Data/p610.txt", header = TRUE)
# 
# ############### Los ahorcados
# 
# X1 = as.data.frame(read.csv("data/LA2.csv",header = T))
# # y <- p610[,5:9]                           # Grab counts
# y <-  as.data.frame(X1[1:nrow(X1),2:11])
# y[y > 1] <- 1                             # Counts to det-nondetections
# C <- sum(apply(y, 1, max)) ; print(C)     # Number of observed species
# table(apply(y, 1, sum))                   # Capture-frequencies
# 
# 
# 
# # Use Augmented data set from Xaug
# # y<- Xaug
# # Augment data set by 150 potential individuals
# nz <- 150
# aug<-array(0,c(nz, ncol(y)))
# colnames(aug)<-colnames(y)
# yaug.la <- rbind(y, aug)
# 
# # Specify model in BUGS language
# sink("M_tbh.jags")
# cat("
#     model {
#     
#     # Priors
#     omega ~ dunif(0, 1)
#     for  (j in 1:T){
#     alpha[j] <- log(mean.p[j] / (1-mean.p[j])) # Define logit 
#     mean.p[j] ~ dunif(0, 1)   # Detection intercepts
#     }
#     gamma ~ dnorm(0, 0.01)
#     tau <- 1 / (sd * sd)
#     sd ~ dunif(0, 3)
#     
#     # Likelihood
#     for (i in 1:M){
#     z[i] ~ dbern(omega)
#     eps[i] ~ dnorm(0, tau)T(-16, 16)
#     
#     # First occasion: no term for recapture (gamma)
#     y[i,1] ~ dbern(p.eff[i,1])
#     p.eff[i,1] <- z[i] * p[i,1]
#     p[i,1] <- 1 / (1 + exp(-lp[i,1]))
#     lp[i,1] <- alpha[1] + eps[i]
#     
#     # All subsequent occasions: includes recapture term (gamma)
#     for (j in 2:T){
#     y[i,j] ~ dbern(p.eff[i,j])
#     p.eff[i,j] <- z[i] * p[i,j]
#     p[i,j] <- 1 / (1 + exp(-lp[i,j]))   
#     lp[i,j] <- alpha[j] + eps[i] + gamma * y[i,(j-1)]
#     } #j
#     } #i
#     
#     # Derived quantities
#     N <- sum(z[])
#     } 
#     ",fill = TRUE)
# sink()
# 
# # Bundle data
# win.data.la <- list(y = as.matrix(yaug.la), M = nrow(yaug.la), T = ncol(yaug.la))
# 
# # Initial values
# inits <- function() list(z = rep(1, nrow(yaug.la)), sd = runif(1, 0.1, 0.9))
# 
# # Parameters monitored
# params <- c("N", "mean.p", "gamma", "sd", "omega", "alpha")
# 
# # MCMC settings
# ni <- 50000
# nt <- 4
# nb <- 10000
# nc <- 3
# 
# # Call JAGS from R (BRT 24 min)
# out.la <- jags(win.data.la, inits, params, "M_tbh.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
# 
# # Summarize posteriors and plot posterior for N
# print(out.la, dig = 3)
# par(mfrow = c(2,1))
# hist(out.la$BUGSoutput$sims.list$N, breaks = 25, col = "gray", main = "", xlab = "Community size LA", las = 1, xlim = c(10, 35), freq = FALSE)
# abline(v = C, col = "red", lwd = 4)
# 
# 
# ########################## Otro sitio
# 
# X2 = as.data.frame(read.csv("data/PS2.csv",header = T))
# # y <- p610[,5:9]                           # Grab counts
# y <-  as.data.frame(X2[1:nrow(X2),2:10])
# y[y > 1] <- 1                             # Counts to det-nondetections
# C <- sum(apply(y, 1, max)) ; print(C)     # Number of observed species
# table(apply(y, 1, sum))                   # Capture-frequencies
# 
# 
# 
# # Use Augmented data set from Xaug
# # y<- Xaug
# # Augment data set by 150 potential individuals
# nz <- 150
# aug<-array(0,c(nz, ncol(y)))
# colnames(aug)<-colnames(y)
# yaug.ps <- rbind(y, aug)
# 
# 
# # Bundle data
# win.data.ps <- list(y = as.matrix(yaug.ps), M = nrow(yaug.ps), T = ncol(yaug.ps))
# 
# # Initial values
# inits <- function() list(z = rep(1, nrow(yaug.ps)), sd = runif(1, 0.1, 0.9))
# 
# # Call JAGS from R (BRT 24 min)
# out.ps <- jags(win.data.ps, inits, params, "M_tbh.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
# 
# # Summarize posteriors and plot posterior for N
# print(out.ps, dig = 3)
# #par(mfrow = c(2,1))
# hist(out.ps$BUGSoutput$sims.list$N, breaks = 20, col = "gray", main = "", xlab = "Community size PS", las = 1, xlim = c(10, 35), freq = FALSE)
# abline(v = C, col = "red", lwd = 4)
# 
# 
# #######################
# ### END 
# #######################
# 
# #### Graficar acumulacion
# 
# alpha.post.ps = out.ps$BUGSoutput$sims.matrix[,"alpha"]
# sigmaU.post.ps = out.ps$BUGSoutput$sims.matrix[,"gamma"]
# N.post.ps = out.ps$BUGSoutput$sims.matrix[,"N"]
# 
# nsites = 30
# CumNumSpeciesPresent(nsites, alpha.post.ps, sigmaU.post.ps, N.post.ps)
# 
# #### Graficar acumulacion
# 
# alpha.post.la = out.la$BUGSoutput$sims.matrix[,"alpha"]
# sigmaU.post.la = out.la$BUGSoutput$sims.matrix[,"gamma"]
# N.post.la = out.la$BUGSoutput$sims.matrix[,"N"]
# 
# nsites = 30
# CumNumSpeciesPresent(nsites, alpha.post.la, sigmaU.post.la, N.post.la)
# 
# 
# 
# 
# # Define model
# sink("M0.jags")
# cat("
#     model {
#     
#     # Priors
#     omega ~ dunif(0, 1)
#     p ~ dunif(0, 1)
#     
#     # Likelihood
#     for (i in 1:M){
#     z[i] ~ dbern(omega)
#     for (j in 1:T){
#     y[i,j] ~ dbern(p.eff[i,j])
#     p.eff[i,j] <- z[i] * p
#     } #j
#     } #i
#     
#     # Derived quantities
#     N <- sum(z[])
#     } # end model
#     ",fill = TRUE)
# sink()
# 
# # Initial values
# inits <- function() list(z = rep(1, nrow(yaug)))
# 
# # Define parameters to be monitored
# params <- c("N", "p", "omega")
# 
# # MCMC settings
# ni <- 5000#0
# nt <- 4
# nb <- 1000#0
# nc <- 3
# 
# # Call JAGS from R (BRT 1 min)
# out0 <- jags(win.data, inits, params, "M0.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,  working.directory = getwd())
# 
# 
# # Summarize posteriors and plot posterior for N
# print(out0, dig = 3)
# par(mfrow = c(1,2))
# hist(out0$BUGSoutput$sims.list$N, breaks = 10, col = "gray", main = "", xlab = "Community size", las = 1, xlim = c(15, 40), freq = FALSE)
# abline(v = C, col = "black", lwd = 3)
# 
# 
