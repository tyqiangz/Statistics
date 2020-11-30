###### Tutorial 8 codes

### Problem 1
## (a) SIR
# first rerun the old codes from Tutorial 6
require(mvtnorm)
s <- 8
n <- 15
t <- 15962989
t1 <- 237217
df <- 4 

logpost <- function(theta, s, n, t, t1){
	theta1 <- theta[1]
	theta2 <- theta[2]
	return((1-s)*theta1 + theta2 - (t-n*t1)*exp(-theta1) - n*exp(theta2-theta1))
}
(out <- optim(par=c(14.5,11), fn=logpost, hessian=TRUE,
              control=list(fnscale=-1), s=s, n=n, t=t, t1=t1))		
(post.mode <- out$par)
(post.cov <- -solve(out$hessian))


# define the difference function
diff <- function(theta, s, n, t, t1, post.mode, post.cov, df){
	logpost(theta,s,n,t,t1) - dmvt(theta,delta=post.mode,sigma=2*post.cov,df=df)
}

set.seed(4234)	
S <- 10^4
# importance sampling		
theta.samples <- rmvt(S,delta=post.mode,sigma=2*post.cov,df=df)
logw <- numeric(S)
for (i in 1:S) {
	logw[i] <- diff(theta.samples[i,], s, n, t, t1, post.mode, post.cov, df)
}
w <- exp(logw - max(logw))
W <- w/sum(w)
hist(W, breaks=20)		# histogram of weights; check for anomalies
# SIR
indices <- sample(1:S, size=S, prob=W, replace=TRUE)
theta.SIR <- theta.samples[indices,]

# contour plot with SIR samples overlayed
theta1.grid <- seq(from=13, to=20, by=.1)
theta2.grid <- seq(from=-10, to=18, by=.1)
theta.grid <- expand.grid(theta1.grid, theta2.grid)			  
grid.logpost <- apply(theta.grid, 1, logpost, s=s, n=n, t=t, t1=t1) 			  
contour(theta1.grid, theta2.grid, matrix(grid.logpost,nrow=length(theta1.grid),ncol=length(theta2.grid)),
        col="blue", nlevels=800, lwd=2, labcex=0.9, 
		xlab=expression(theta[1]), ylab=expression(theta[2]))
points(theta.SIR[,1], theta.SIR[,2], col=rgb(1,0,0,alpha=0.2),pch=19)		


## (b) 
t0 <- 10^6
Rt0 <- function(theta, t0, t1){
	mu <- t1 - exp(theta[,2])
	beta <- exp(theta[,1])
	R <- exp(-(t0-mu)/beta)
	return(R)
}
R.samples <- Rt0(theta.SIR,t0,t1)
mean(R.samples)
sd(R.samples)



### Problem 2
## (a) Random walk Metropolis
# log posterior function
logpost <- function(eta) {
	125*log(2+(exp(eta)/(1+exp(eta))))-74*log(1+exp(eta))+35*eta	
}

# find the posterior mode using optim() with method="Brent"
(out <- optim(par=0, fn=logpost, hessian=TRUE, control=list(fnscale=-1), 
              method="Brent", lower=-10, upper=10))
(post.mode <- out$par)
(post.var <- -1/out$hessian)

# random walk Metropolis
require(LearnBayes)
T <- 10^4
proposal <- list(var=post.var,scale=3)
set.seed(4234)
fit1 <- rwmetrop(logpost,proposal,start=post.mode,m=T)
fit1$accept

# histogram of RW Metropolis draws and normal approximation
# normal curve is overlayed in the histogram
eta.grid <- seq(from=min(fit1$par), to=max(fit1$par), length=1000)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
hist(fit1$par, freq=F, breaks=20, col="cyan", 
     main="RW Metropolis Draws", xlab=expression(eta), 
	 ylab=expression(paste("p(",eta,"|y)")))
lines(eta.grid, dnorm(eta.grid,mean=post.mode,sd=sqrt(post.var)), 
      lty=2,lwd=2,col="red")

# 95% Bayesian CI of theta
CI <- quantile(fit1$par[5001:T], probs=c(0.025,0.975))
exp(CI)/(1+exp(CI))


## (b) Independence Metropolis
proposal2 <- list(mu=post.mode,var=9*post.var)
set.seed(4234)
fit2 <- indepmetrop(logpost,proposal2,start=post.mode,m=T)
fit2$accept

# histogram of RW Metropolis draws and normal approximation
# normal curve is overlayed in the histogram
eta.grid <- seq(from=min(fit2$par), to=max(fit2$par), length=1000)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
hist(fit2$par, freq=F, breaks=20, col="cyan", 
     main="Independence Metropolis Draws", xlab=expression(eta), 
	 ylab=expression(paste("p(",eta,"|y)")))
lines(eta.grid, dnorm(eta.grid,mean=post.mode,sd=sqrt(post.var)), 
      lty=2,lwd=2,col="red")

# 95% Bayesian CI of theta
CI <- quantile(fit2$par[5001:T], probs=c(0.025,0.975))
exp(CI)/(1+exp(CI))



### Problem 3
## (b)
y <- c(15,11,14,17,5,11,10,4,8,10,7,9,11,3,6,1,1,4)
n <- length(y)
x <- 1:n

# log posterior function
logpost <- function(theta,y){
	n <- length(y)
	beta0 <- theta[1]
	beta1 <- theta[2]
	lp <- beta0 + beta1*(1:n)
	L <- sum(y*lp - exp(lp))
	return(L)
}

# first fit the generalized linear model (Poisson regression)
fit.glm <- glm(y~x,family=poisson)
summary(fit.glm)
theta.mle <- fit.glm$coefficients	# MLE of theta used as initial values


(out <- optim(par=theta.mle,fn=logpost,hessian=TRUE,
              control=list(fnscale=-1),y=y))
(post.mode <- out$par)
(post.cov <- -solve(out$hessian))


## (c)
# (i) Random walk Metropolis
require(LearnBayes)
require(coda)
T <- 10^4
proposal1 <- list(var=post.cov,scale=2)
start <- post.mode
set.seed(4234)
fit1 <- rwmetrop(logpost,proposal1,start,T,y)
fit1$accept


colnames(fit1$par) <- c("beta0","beta1")
plot(mcmc(fit1$par))
autocorr.plot(mcmc(fit1$par))
mean(fit1$par[,2])
sd(fit1$par[,2])

# (ii) Independence Metropolis
proposal2 <- list(mu=post.mode, var=4*post.cov)
set.seed(4234)
fit2 <- indepmetrop(logpost,proposal2,start,T,y)
fit2$accept

colnames(fit2$par) <- c("beta0","beta1")
plot(mcmc(fit2$par))
autocorr.plot(mcmc(fit2$par))
mean(fit2$par[,2])
sd(fit2$par[,2])


## (d) 5%, 50%, 95% quantiles of three methods

# normal approximation
post.sd <- sqrt(diag(post.cov))
qnorm(c(0.05,0.5,0.95),mean=post.mode[1],sd=post.sd[1])
qnorm(c(0.05,0.5,0.95),mean=post.mode[2],sd=post.sd[2])

# RW Metropolis
apply(fit1$par,2,quantile,probs=c(0.05,0.5,0.95))

# Independence Metropolis
apply(fit2$par,2,quantile,probs=c(0.05,0.5,0.95))


