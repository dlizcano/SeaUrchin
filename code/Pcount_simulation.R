


# Simulate data
set.seed(35)
nSites <- 16
nVisits <- 4
x <- rnorm(nSites)               # a covariate
beta0 <- 0
beta1 <- 1
lambda <- exp(beta0 + beta1*x)   # expected counts at each site
N <- rpois(nSites, lambda)       # latent abundance
y <- matrix(NA, nSites, nVisits)
p <- c(0.3, 0.6, 0.8, 0.5)            # detection prob for each visit
for(j in 1:nVisits) {
  y[,j] <- rbinom(nSites, N, p[j])
}

# Organize data
visitMat <- matrix(as.character(1:nVisits), nSites, nVisits, byrow=TRUE)

umf <- unmarkedFramePCount(y=y, siteCovs=data.frame(x=x),
                           obsCovs=list(visit=visitMat))
summary(umf)

# Fit a model
fm1 <- pcount(~visit-1 ~ x, umf, K=50)
fm1

plogis(coef(fm1, type="det")) # Should be close to p


# Empirical Bayes estimation of random effects
(fm1re <- ranef(fm1))
plot(fm1re, subset=site \%in\% 1:25, xlim=c(-1,40))
sum(bup(fm1re))         # Estimated population size
sum(N)                  # Actual population size
