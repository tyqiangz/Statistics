###### Chapter 9 Hierarchical Models ######
### All following codes are for heart transplant mortality data

## read heart transplant data
require(LearnBayes)
data(hearttransplants)
y <- hearttransplants$y
e <- hearttransplants$e

# plot on page 9
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(log(e),y/e, pch=as.character(y),col="blue")

## Equal rates model
a <- sum(y)
b <- sum(e)
i <- 94

# plot on page 18
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
grid <- seq(from=0,to=30,by=1)
plot(dnbinom(grid,size=a,prob=b/(b+e[i])),type="h",
     col="blue",lwd=6,ylab="Density")
abline(v=y[i],col="red",lwd=2)


prob1 <- pnbinom(y,size=a,prob=b/(b+e))
prob2 <- 1-prob1+ dnbinom(y,size=a,prob=b/(b+e))
prob <- apply(cbind(prob1,prob2),1,min)

# plot on page 20
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(log(e),prob)
sum(prob<0.1)

## Exchangeable model
# marginal posterior of alpha and mu

z0 <- 0.53
data <- list(z0=z0,y=y,e=e)

logpost <- function(theta, data){
	y <- data$y
	e <- data$e
	z0 <- data$z0
	n <- length(y)
	alpha <- exp(theta[1])
	mu <- exp(theta[2])
	L <- (-n*lgamma(alpha) + theta[1] -2*log(alpha+z0) 
         + 94*alpha*(theta[1]-theta[2]) + sum(lgamma(alpha+y))
         - sum((alpha+y)*log(alpha/mu+e)))
	return(L)
}

# contour plot on page 29
N <- 100
theta1 <- seq(from=-2, to=16, length.out=N)
theta2 <- seq(from=-8, to=-6, length.out=N)
z <- matrix(0,N,N)
for (i in 1:N){
	for (j in 1:N){
		theta <- c(theta1[i], theta2[j])
		z[i,j] <- logpost(theta,data)    
}}

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
contour(x=theta1,y=theta2,z,col="blue",nlevels=10,
        xlab=expression(theta[1]),ylab=expression(theta[2]),
        cex.axis=1.2,cex.lab=1.4,labcex=1)


# find normal approximation
start <- c(2,-7)
out <- optim(par=start,fn=logpost,hessian=TRUE,
             control=list(fnscale=-1),data=data)
(post.mode <- out$par)
(post.cov <- -solve(out$hessian))
post.sd <- sqrt(diag(post.cov))


# Metropolis within Gibbs 
# (we can also use other algorithms, such as rwmetrop, indepmetrop, etc.)
T <- 10000
set.seed(1)
fit <- gibbs(logpost, start=start,m=T,scale=2.5*post.sd,data=data)
fit$accept

# MCMC output analysis
require(coda)
mcmcobj <- mcmc(fit$par)
colnames(mcmcobj) <- c("theta1", "theta2")

# trace plots and density plots
par(mar=c(3.5,2,1,1))
par(mgp=c(2.1,0.8,0))
plot(mcmcobj)

# autocorrelation plots
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(2,1))
autocorr.plot(mcmcobj,auto.layout=FALSE)

# contour plot with posterior draws overlayed on page 33
N <- 100
theta1 <- seq(from=0, to=8, length.out=N)
theta2 <- seq(from=-7.3, to=-6.6, length.out=N)
z <- matrix(0,N,N)
for (i in 1:N){
	for (j in 1:N){
		theta <- c(theta1[i], theta2[j])
		z[i,j] <- logpost(theta,data)    
}}

par(mfrow=c(1,1))
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
contour(x=theta1,y=theta2,z,col="blue",nlevels=20,
        xlab=expression(theta[1]),ylab=expression(theta[2]),
        cex.axis=1.2,cex.lab=1.5,labcex=0.9)
points(fit$par, col="red")


# simulate rates (lambda)
alpha_draws <- exp(fit$par[,1])
mu_draws <- exp(fit$par[,2]) 
lambda_draws <- matrix(0,T,94)
lambdaCI <- matrix(0,94,2)
for (i in 1:94){
	lambda_draws[,i] <- rgamma(10000,alpha_draws+y[i], alpha_draws/mu_draws+e[i])
	lambdaCI[i,] <- quantile(lambda_draws[,i],c(0.05,0.95))
}

# plot lambda on page 36
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(log(e),y/e, pch=as.character(y),col="blue")
for (i in 1:94){
	lines(rep(log(e[i]),2),lambdaCI[i,])
}


## Shrinkages: plot on page 40
B <- rep(0,94)
for (i in 1:94){
	B[i] <- mean(alpha_draws/(alpha_draws + e[i]*mu_draws))
}
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(log(e),B,ylab="shrinkage")


# comparing hospitals 
lambda_post.mean <- apply(lambda_draws,2,mean)
order(lambda_post.mean)[1]

# all possible pairs
choose(94,2)
compare <- matrix(0,94,94)
for (i in 1:93){
	for (j in (i+1):94){
		compare[i,j] <- sum(lambda_draws[,i] < lambda_draws[,j])/T
		compare[j,i] <- 1 - compare[i,j]
}}
compare[1:24,85]


## Sensitivity analysis

# SIR algorithm: get samples from the new posterior,
# by applying SIR to the draws from the old posterior
theta1_draws <- fit$par[,1]
logprior <- function(theta1,z0){
	log(z0) + theta1 - 2*log(z0 + exp(theta1))
}

logw <- logprior(theta1_draws,5) - logprior(theta1_draws,0.53)
w <- exp(logw-max(logw))
W <- w/sum(w)
idx <- sample(1:T, size=10000,replace=TRUE, prob=W)
theta1_new <- theta1_draws[idx]

# comparison of two posterior densities, plot on page 48
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
theta1_grid <- seq(from=-3,to=5,by=0.1)
plot(theta1_grid,exp(logprior(theta1_grid,0.53)),
     type="l",lty=2,xlab=expression(theta[1]),ylab="density",
     ylim=c(0,0.82),lwd=2)
points(theta1_grid,exp(logprior(theta1_grid,5)),
       type="l",lty=2,col="red",lwd=2)
points(density(theta1_draws),type="l",lwd=2)
points(density(theta1_new),type="l",lwd=2,col="red")
legend("topleft",legend=c("prior(z0=0.53)","prior(z0=5)",
       "posterior(z0=0.53)","posterior(z0=5)"),lty=c(2,2,1,1),
       col=c("black","red","black","red"),lwd=2)


## Posterior predictive model checking
set.seed(1)
y94_draws <- rpois(T, e[94]*lambda_draws[,94])

# histogram on page 51
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
hist(y94_draws)
abline(v=y[94],col="red",lwd=2)


# probabilities of extremes
prob.exchange <- rep(0,94)
for (i in 1:94){
	yi_draws <- rpois(T, e[i]*lambda_draws[,i])
	prob1 <- sum(yi_draws <= y[i])/T
	prob2 <- sum(yi_draws >= y[i])/T
	prob.exchange[i] <- min(prob1,prob2)
}

# plot on page 53
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(prob, prob.exchange, xlab="P(extreme),equal means",
  ylab="P(extreme),exchangeable")
abline(0,1)

sum(prob.exchange<0.1)
	   
	   
	   
	   
