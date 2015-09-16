# Fit the abundance-based community model in JAGS
# The data are a subset of the removal count data from
# Chandler et al. (Conservation Biology, In press)


load("data-comm-abund.RData")
library(rjags)
ls()




# The model
sink("model-communityN.txt")
cat("
model {
int.mu.beta0 ~ dnorm(0, 0.01)
int.mu.beta1 ~ dnorm(0, 0.01)
int.mu.beta2 ~ dnorm(0, 0.01)
int.mu.beta3 ~ dnorm(0, 0.01)
int.mu.beta4 ~ dnorm(0, 0.01)
x.mu.beta0 ~ dnorm(0, 0.01)
x.mu.beta1 ~ dnorm(0, 0.01)
x.mu.beta2 ~ dnorm(0, 0.01)
x.mu.beta3 ~ dnorm(0, 0.01)
x.mu.beta4 ~ dnorm(0, 0.01)
sig.beta0 ~ dunif(0, 5)
tau.beta0 <- 1/(sig.beta0*sig.beta0)
sig.beta1 ~ dunif(0, 5)
tau.beta1 <- 1/(sig.beta1*sig.beta1)
sig.beta2 ~ dunif(0, 5)
tau.beta2 <- 1/(sig.beta2*sig.beta2)
sig.beta3 ~ dunif(0, 5)
tau.beta3 <- 1/(sig.beta3*sig.beta3)
#mu.beta4 ~ dnorm(0, 0.01)
sig.beta4 ~ dunif(0, 5)
tau.beta4 <- 1/(sig.beta4*sig.beta4)
mu.logitp ~ dnorm(0, 0.1)
sig.logitp ~ dunif(0, 5)
tau.logitp <- 1/(sig.logitp*sig.logitp)
alpha1 ~ dnorm(0, 0.01)
alpha2 ~ dnorm(0, 0.01)
alpha3 ~ dnorm(0, 0.01)
alpha4 ~ dnorm(0, 0.01)
alpha5 ~ dnorm(0, 0.01)
psi ~ dbeta(1,1) # For data-augmentation (we called this omega)
for(g in 1:5) { # 5 forest-dependence classes. Guild-level
  mu.beta0[g] <- int.mu.beta0 + x.mu.beta0*(g-1)
  mu.beta1[g] <- int.mu.beta1 + x.mu.beta1*(g-1)
  mu.beta2[g] <- int.mu.beta2 + x.mu.beta2*(g-1)
  mu.beta3[g] <- int.mu.beta3 + x.mu.beta3*(g-1)
  mu.beta4[g] <- int.mu.beta4 + x.mu.beta4*(g-1)
  }
for(i in 1:(nspp+nz)) {
  w[i] ~ dbern(psi)
  # Species-specific coefficients
  beta0[i] ~ dnorm(mu.beta0[depindex[i]], tau.beta0)
  beta1[i] ~ dnorm(mu.beta1[depindex[i]], tau.beta1)
  beta2[i] ~ dnorm(mu.beta2[depindex[i]], tau.beta2)
  beta3[i] ~ dnorm(mu.beta3[depindex[i]], tau.beta3)
  beta4[i] ~ dnorm(mu.beta4[depindex[i]], tau.beta4)
  alpha0[i] ~ dnorm(mu.logitp, tau.logitp)
  for(j in 1:32) {
     log(lambda[i,j]) <- beta0[i] +
                         beta1[i]*ioc[j] +
                         beta2[i]*shade[j] +
                         beta3[i]*secondary[j] +
                         beta4[i]*distDivide[j]
     logit(p[i,j,1]) <- alpha0[i] + alpha1*nethours[j,1] +
         alpha2*wind[j,1] + alpha3*ioc[j] + alpha4*shade[j] +
         alpha5*secondary[j]
     logit(p[i,j,2]) <- alpha0[i] + alpha1*nethours[j,2] +
         alpha2*wind[j,2] + alpha3*ioc[j] + alpha4*shade[j] +
         alpha5*secondary[j]
     logit(p[i,j,3]) <- alpha0[i] + alpha1*nethours[j,3] +
         alpha2*wind[j,3] + alpha3*ioc[j] + alpha4*shade[j] +
         alpha5*secondary[j]
     N[i,j] ~ dpois(lambda[i, j] * w[i])
     y[i,j,1] ~ dbin(p[i,j,1], N[i,j])
     N2[i,j] <- N[i,j]-y[i,j,1]
     y[i,j,2] ~ dbin(p[i,j,2], N2[i,j])
     N3[i,j] <- N2[i,j]-y[i,j,2]
     y[i,j,3] ~ dbin(p[i,j,3], N3[i,j])
     }
  # Total abundance in each habitat type
  Nh[i,1] <- N[i,] %*% shade
  Nh[i,2] <- N[i,] %*% ioc
  Nh[i,3] <- N[i,] %*% secondary
  Nh[i,4] <- N[i,] %*% primary
  }
  S <- sum(w[])
}
", fill=TRUE)
sink()








Ni <- array(0L, c(nrow(jdat4$y), 32))
Ni[1:75,] <- apply(jdat4$y[1:75,,], c(1,2), sum)+2
str(Ni)

jinit <- function() list(N=Ni,
                         w=as.integer( apply(Ni>0, 1, any)),
#                         w=rep(1L, nrow(Ni)),
                         int.mu.beta0=rnorm(1), x.mu.beta0=rnorm(1),
                         int.mu.beta1=rnorm(1), x.mu.beta1=rnorm(1),
                         int.mu.beta2=rnorm(1), x.mu.beta2=rnorm(1),
                         int.mu.beta3=rnorm(1), x.mu.beta3=rnorm(1),
                         int.mu.beta4=rnorm(1), x.mu.beta4=rnorm(1),
                         sig.beta0=runif(1), sig.beta1=runif(1),
                         sig.beta2=runif(1), sig.beta3=runif(1),
                         sig.beta4=runif(1),
                         mu.logitp=rnorm(1), sig.logitp=runif(1, 1, 2),
                         alpha1=rnorm(1), alpha2=rnorm(1))
str(jinit())


jpars <- c(paste("int.mu.beta", 0:4, sep=""),
           paste("x.mu.beta", 0:4, sep=""),
           paste("sig.beta", 0:4, sep=""),
           "mu.logitp", "sig.logitp",
           "alpha1", "alpha2", "alpha3", "alpha4", "alpha5",
           "S")




jm <- jags.model("model-communityN.txt", jdat4, jinit, n.adapt=500,
                 n.chains=1)
jc1 <- coda.samples(jm, jpars, n.iter=1000, thin=1)


plot(jc1, ask=TRUE)
summary(jc1, start=1001)


jm4state <- jm4$state()
str(jm4state)


#save(jdat4, jm4state,

save.image("err.RData")



jc2 <- coda.samples(jm, jpars, n.iter=5000, thin=1)

plot(jc2, ask=TRUE)
summary(jc2, start=1001)



jc3 <- coda.samples(jm, jpars, n.iter=5000, thin=1)

plot(jc3, ask=TRUE)
summary(jc3, start=1001)




sort(sapply(ls(), function(x) object.size(get(x))))






# Monitor the population size parameters
jcNh <- coda.samples(jm, "Nh", n.iter=100, thin=1)






# Compute posteriors of community parameters from abundance matrix
cstats <- function(mc) {
    if(!is.mcmc.list(mc))
        stop("mc should have class 'mcmc.list'")
    nspp <- 75 # Number of observed species. Change this for your dataset!
    forspp. <- as.integer(jdat4$depindex < 3)
    nc <- nchain(mc)
    ni <- niter(mc)*nc
    m <- as.matrix(mc)
    N <- array(m, c(ni, nspp+jdat4$nz, 4)) # abundance
    dimnames(N) <- list(NULL, NULL, #forestDepend$Scientific,
                        c("shade", "ioc", "sf", "pf"))
    R <- apply(N, c(1,3), function(x) x/sum(x)) # Relative abundance
    R <- aperm(R, c(2,1,3))
    Rno0 <- R
    Rno0[R==0] <- .Machine$double.xmin # avoid log(0)
    Z <- ifelse(N>0, 1, 0)
    S <- apply(Z, c(1,3), sum)
    # A=shared, B=in x but not PF, C=in PF but not x
    A <- t(apply(Z, 1, function(x)
                 (x[,1:3]==1) & (x[,4]==1)))
    A <- array(A, c(ni, nspp+jdat4$nz, 3)) # double check order
    A[] <- as.integer(A)
    B <- t(apply(Z, 1, function(x)
                 (x[,1:3]==1) & (x[,4]==0)))
    B <- array(B, c(ni, nspp+jdat4$nz, 3)) # double check order
    B[] <- as.integer(B)
    C <- t(apply(Z, 1, function(x)
                 (x[,1:3]==0) & (x[,4]==1)))
    C <- array(C, c(ni, nspp+jdat4$nz, 3)) # double check order
    C[] <- as.integer(C)
    # D is the union of x and PF
    D <- t(apply(Z, 1, function(x)
                 (x[,1:3]==1) | (x[,4]==1)))
    D <- array(D, c(ni, nspp+jdat4$nz, 3))
    D[] <- as.integer(D)
    S2 <- apply(Z, c(1,3), function(x)
                x %*% forspp.)
    H <- apply(Rno0, c(1,3), function(x) -sum(x * log(x)))
    H2 <- apply(Rno0, c(1,3), function(x)
                -((x * log(x)) %*% forspp.))
    J <- J2 <- Jm <- CJ <- CJ2 <- CJm <- Uniq <- Uniq2 <- UniqM <-
        Share <- Share2 <- ShareM <-
        matrix(NA, ni, 3, dimnames=list(NULL, c("shade", "ioc", "sf")))
    for(i in 1:ni) {
        Union <- colSums(D[i,,])
        Union2 <- t(D[i,,] ) %*% forspp.
        Uniq[i,] <- colSums(C[i,,]) #/ Union
        Uniq2[i,] <- (t(C[i,,]) %*% forspp.) #/ Union2
        Share[i,] <- colSums(A[i,,])
        Share2[i,] <- t(A[i,,]) %*% forspp.
        J[i,] <- colSums(A[i,,]) / Union
        J2[i,] <- (t(A[i,,]) %*% forspp.) / Union2
        U <- colSums(R[i,,1:3]*A[i,,]) # Relative abund for shared spp in X
        V <- colSums(R[i,,4]*A[i,,])  # Relative abund for shared spp in PF
        CJ[i,] <- (U*V) / (U + V + (U*V))
        U2 <- colSums(R[i,,1:3]*A[i,,]*forspp.)
        V2 <- colSums(R[i,,4]*A[i,,]*forspp.)
        CJ2[i,] <- (U2*V2) / (U2 + V2 + (U2*V2))
    }
    # Order columns as Shade, IOC, SF, PF
    return(list(S=S,    # Species richness of all species
                S2=S2,  # Species richness of forest-dependents
                H=H,    # Shannon diversity of all spp
                H2=H2,
                J=J,    # Jaccard similarity index
                J2=J2,
                CJ=CJ,  # Chao-Jaccard index
                CJ2=CJ2,
                Uniq=Uniq, # Unique to PF
                Uniq2=Uniq2,
                Share=Share, # Shared species
                Share2=Share2)) #,
}



v <- cstats(jcNh)
str(v)



boxplot(v$S)
boxplot(v$S2)

boxplot(v$H)
boxplot(v$H2)

boxplot(v$CJ)
boxplot(v$CJ2)

boxplot(v$J)
boxplot(v$J2)

boxplot(v$Uniq)
boxplot(v$Uniq2)

boxplot(v$Share)
boxplot(v$Share2)




















