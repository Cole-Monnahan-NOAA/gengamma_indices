### Demo script of applying rescaling bootstrap estimators from
### Smith (1997) CJFAS.

### Some functions for fittin the GGD in R
#' Convert a mean to a 'mu' parameter
#' @param x An observed value
#' @param mean The mean of the distribution
#' @param sigma,Q
mean2mu <- function(mean,sigma,Q){
  if(Q==0) return(log(mean)) ## lognormal case
  k <- 1/(Q*Q)
  beta <- Q/sigma;
  log_theta <- log(mean) - lgamma( (k*beta+1)/beta ) + lgamma( k );
  mu <- log_theta + log(k) / beta
  return(mu)
}
dgengamma_mean <- function(x, mean, sigma, Q, log=FALSE)
  flexsurv::dgengamma(x=x, mu=mean2mu(mean, sigma,Q),
                      sigma=sigma, Q=Q, log=log)
nll_ggd <- function(pars, samples){
  -sum(dgengamma_mean(x=samples, exp(pars[1]), sigma=exp(pars[2]), Q=pars[3],
                 log=TRUE))
}
fit_ggd <- function(pars, samples){
  opt <- optim(par=c(3,1,.1), fn=nll_ggd, samples=samples)
  return(c(exp(opt$par[1]), exp(opt$par[2]), opt$par[3]))
}
fit_lognormal <- function(pars, samples){
  y <- samples
  nll_ln <- function(pars, y) -sum(dlnorm(x=y, meanlog=pars[1], sdlog=exp(pars[2]), log=TRUE))
  opt <- optim(par=pars, fn=nll_ln, y=y)
  return(c(opt$par[1], exp(opt$par[2])))
}



## global bootstrapping settings
nboot <- 5000                  # no. of boot replicates
L <- 20                        # no. of strata
haul_scalar <- 5               # scale up how many possible hauls
sim.seed <- 1541               # for simulating data
boot.seed <- 3421              # for bootstrapping
moffset <- 3                   # how many fewer samples for rescaled?
outlier <- FALSE                # include a giant outlier tow?


## simulate some data
set.seed(sim.seed)
Nh <- haul_scalar*sample( 10:50, size=L, replace=TRUE) # no. hauls per stratum
nh <- ceiling(Nh*.1)                         # observed no. hauls
N <- sum(Nh)                               # total hauls
hauls <- lapply(nh, function(i) rlnorm(i, 4,.5))
## simulate a huge unexpected tow?
if(outlier) hauls[[3]][1] <- max(hauls[[3]])*15
W <- runif(L, .2,.6)                      # area weightings?
W <- W/sum(W)
xseq <- seq(0,200, len=1000)
ybar <- sapply(hauls, mean)
varbar <- sapply(hauls, function(x) var(x))
par(mfrow=c(4,5))
## for(i in 1:L){
##   hist(hauls[[i]], freq=FALSE)
##   lines(xseq, dnorm(xseq, ybar[i], sqrt(varbar[i])))
## }
## DB estimates of total biomass of original data set
yhat <- sum(W*ybar)                     # eqn (1)
varhat <- sum(Nh/N^2*(Nh-nh)*varbar/nh) # eqn (2)



## Naieve bootstrapping
set.seed(boot.seed)
varbarboots <- yhatboots_naive <- rep(NA,nboot)
for(i in 1:nboot){
  haulsboot <- lapply(hauls, function(x) sample(x, size=length(x), replace=TRUE))
  ybootbar <- sapply(haulsboot, mean)
  yhatboots_naive[i] <- sum(W*ybootbar)
}
## replicate estimates of total biomass and variance from boots
yhat_naive <- mean(yhatboots_naive)             # eqn 3
varhat_naive <- sum((yhatboots_naive-yhat_naive)^2)/(nboot-1) # eqn 4
ggd_pars_naive <- fit_ggd(pars=c(log(yhat_naive), log(varhat_naive), -.1), samples=yhatboots_naive)
ln_pars_naive <- fit_lognormal(pars=c(log(yhat_naive), .1), samples=yhatboots_naive)

## scaled bootstrapping
set.seed(boot.seed)
                                        # to take
if(any(nh<=moffset)) stop("too few hauls for moffset=", moffset)
fh <- nh/Nh
mh <- nh-moffset
varbarboots <- yhatboots_rescaled<- rep(NA,nboot)
scale <- sqrt(mh*(1-fh)/(nh-1))         # scaling term by stratum
haulboot <- haulsbootscaled <- list()
for(i in 1:nboot){
  ## loop over each stratum
  for(j in 1:L){
    ## sample smaller size (mh)
    haulsboot[[j]] <- sample(hauls[[j]], size=mh[j], replace=TRUE)
    ## now "scale" the samples, where ybar is the original mean
    haulsbootscaled[[j]] <- ybar[j] + scale[j]*(haulsboot[[j]]-ybar[j])
  }
  ybootbar <- sapply(haulsbootscaled, mean)
  yhatboots_rescaled[i] <- sum(W*ybootbar)
}
yhat_rescaled<- mean(yhatboots_rescaled)
varhat_rescaled<- sum((yhatboots_rescaled-yhat_rescaled)^2)/(nboot-1)
ggd_pars_rescaled <- fit_ggd(pars=c(log(yhat_rescaled), log(varhat_rescaled), -.1),
                         samples=yhatboots_rescaled)
ln_pars_rescaled <- fit_lognormal(pars=c(log(yhat_rescaled), .1), samples=yhatboots_rescaled)





par(mfrow=c(3,1), mgp=c(3,.5,0), tck=-.01, mar=c(3,3,2,.5))
breaks <- 40
xseq <- seq(min(yhatboots_naive)/1.005, max(yhatboots_naive)*1.005, len=1000)
hist(yhatboots_naive, freq=FALSE, main='Naive bootstrap', breaks=breaks, xlim=range(xseq))
lines(xseq, dnorm(xseq, yhat_naive, sqrt(varhat_naive)), lwd=2)
lines(xseq, dnorm(xseq, yhat, sqrt(varhat)), col=2, lwd=2)
lines(xseq, dgengamma_mean(xseq, ggd_pars_naive[1], ggd_pars_naive[2], ggd_pars_naive[3]), col=4, lwd=2)
lines(xseq, dlnorm(xseq, ln_pars_naive[1], ln_pars_naive[2]), col=3, lwd=2)
legend("right", legend=c("Original DB", "Boot DB", "lognormal MLE", "gengamma MLE"),
       col=1:4, lty=1)
hist(yhatboots_rescaled, freq=FALSE, main='Rescaled bootstrap', breaks=breaks, xlim=range(xseq))
lines(xseq, dnorm(xseq, yhat_rescaled, sqrt(varhat_rescaled)), lwd=2)
lines(xseq, dnorm(xseq, yhat, sqrt(varhat)), col=2, lwd=2)
lines(xseq, dgengamma_mean(xseq, ggd_pars_rescaled[1], ggd_pars_rescaled[2], ggd_pars_rescaled[3]), col=4, lwd=2)
lines(xseq, dlnorm(xseq, ln_pars_rescaled[1], ln_pars_rescaled[2]), col=3, lwd=2)
legend("right", legend=c("Original DB", "Boot DB", "lognormal MLE", "gengamma MLE"),
       col=1:4, lty=1)
plot(xseq,
     dgengamma_mean(xseq, ggd_pars_naive[1], ggd_pars_naive[2], ggd_pars_naive[3]),
     col=1, lwd=2, type='l', ylab='density')
lines(xseq, dgengamma_mean(xseq, ggd_pars_rescaled[1], ggd_pars_rescaled[2], ggd_pars_rescaled[3]), col=2, lwd=2)
legend('topright', legend=c('Naive', 'Rescaled'), col=1:2, lty=1)
