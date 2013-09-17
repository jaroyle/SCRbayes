\name{ch9secrYjags}
\alias{ch9secrYjags}
\title{
  Fit IPP using secr and JAGS
}
\description{
  The analysis of is of simulated data with a single state-space
  covariate. See Chapter 9 for an application using real data
}
\details{
  See Chapter 9 or is it 10?
}
\seealso{
  \code{\link{scrIPP}}
}
\examples{

\dontrun{

library(secr)
library(rjags)

data(ch9simData)

ch <- ch9simData$ch.secr
msk <- ch9simData$spcov.secr


# SECR analysis

secr1 <- secr.fit(ch, model=D~elev, mask=msk)

region.N(secr1, se.N=TRUE)





# JAGS analysis

# JAGS model
sink("ippDiscrete.txt")
cat("
model{
sigma ~ dunif(0, 1)
lam0 ~ dunif(0, 5)
beta ~ dnorm(0,0.1)
psi ~ dbeta(1,1)

for(j in 1:nPix) {
  theta[j] <- exp(beta*elevation[j])
  probs[j] <- theta[j]/sum(theta[])
}

for(i in 1:M) {
  w[i] ~ dbern(psi)
  s[i] ~ dcat(probs[])
  x0g[i] <- Sgrid[s[i],1]
  y0g[i] <- Sgrid[s[i],2]
  for(j in 1:ntraps) {
    dist[i,j] <- sqrt(pow(x0g[i]-grid[j,1],2) + pow(y0g[i]-grid[j,2],2))
    lambda[i,j] <- lam0*exp(-dist[i,j]*dist[i,j]/(2*sigma*sigma)) * w[i]
    y[i,j] ~ dpois(lambda[i,j])
    }
  }

N <- sum(w[])
D <- N/1 # unit square
}

", fill=TRUE)
sink()



modfile <- "ippDiscrete.txt"

jags.data <- with(ch9simData, {
    list(y=ch.jags, elevation=drop(spcov.jags$elev),
            nPix=nrow(spcov.jags),
            M=nrow(ch.jags), ntraps=nrow(traps),
            Sgrid=as.matrix(spcov.jags[,1:2]),
            grid=traps)
    })
str(jags.data)

init <- function() {
    list(sigma=runif(1), lam0=runif(1), beta=rnorm(1),
         s=sample.int(jags.data$nPix, jags.data$M, replace=TRUE),
         w=rep(1, jags.data$M), psi=1)
}
str(init())

pars <- c("sigma", "lam0", "beta", "N")

# Obtain posterior samples. This takes a few minutes
# Compile and adapt
set.seed(03453)
jm <- jags.model(modfile, jags.data, init, n.chains=2, n.adapt=200)
# MCMC
jags1 <- coda.samples(jm, pars, n.iter=2200)

plot(jags1)
summary(window(jags1, start=1001))


unlink(modfile)


}



}
