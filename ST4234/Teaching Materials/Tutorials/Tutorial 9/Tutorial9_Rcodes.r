####### Tutorial 9 codes

require(LearnBayes)
require(coda)

### Problem 1
## (b) Gibbs sampler
y <- c(125,18,20,34)
# y <- c(14,0,1,5)
T <- 10^4
theta.trace <- numeric(T)
z.trace <- numeric(T)
theta.old <- 0.5	# initialize theta

set.seed(4234)
for (i in 1:T) {
	z.new <- rbinom(1,y[1],2/(2+theta.old))
	theta.new <- rbeta(1,y[1]-z.new+y[4]+1,y[2]+y[3]+1)
	theta.trace[i] <- theta.new
	z.trace[i] <- z.new
	theta.old <- theta.new
}

# traceplot and density plot of theta
plot(mcmc(theta.trace[5001:T]))

# posterior mean and sd of theta
mean(theta.trace[5001:T])
sd(theta.trace[5001:T])



### Problem 2
y <- c(9.3,4.9,3.5,26.0,0.6,1.0,3.5,26.9,2.6,20.4,1.0,10.0,1.7,11.3,7.7,
      14.1,24.8,3.8,8.4,1.1,24.5,90.7,16.4,30.7,8.5,5.9,14.7,0.5,99.5,35.2)
n <- length(y)

## (b)
logpost <- function(theta,y){
	lambdaA <- exp(theta[1])
	lambdaB <- exp(theta[2])
	return(sum(log(0.8*exp(-y/lambdaA)/lambdaA + 0.2*exp(-y/lambdaB)/lambdaB)))
}

# contour plot
theta1 <- seq(from=1, to=4, by=0.1)
theta2 <- seq(from=-2, to=8, by=0.1)
z <- matrix(0,length(theta1),length(theta2))
for (i in 1:length(theta1)){
	for (j in 1:length(theta2)){
		theta <- c(theta1[i], theta2[j])
		z[i,j] <- logpost(theta,y)    
}}

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
contour(x=theta1,y=theta2,z,col="blue",nlevels=40,
        xlab=expression(theta[A]),ylab=expression(theta[B]),
        cex.axis=1.1,cex.lab=1.3)

## (c)
(out1 <- optim(par=c(3,0),fn=logpost,hessian=TRUE,control=list(fnscale=-1),y=y))
		
## (d)
(out2 <- optim(par=c(2,4),fn=logpost,hessian=TRUE,control=list(fnscale=-1),y=y))

## (f)
# starting at the 1st mode
post.mode <- out1$par
post.cov <- -solve(out1$hessian)

T <- 10^4
start <- post.mode
proposal <- list(var=post.cov,scale=3)
set.seed(4234)
fit1.1<- rwmetrop(logpost,proposal,start,T,y)
fit1.1$accept

# density plot
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2))
plot(density(fit1.1$par[5001:T,1]),main="",xlab=expression(theta[A]))
plot(density(fit1.1$par[5001:T,2]),main="",xlab=expression(theta[B]))

# starting at the 2nd mode
post.mode <- out2$par
post.cov <- -solve(out2$hessian)

T <- 10^4
start <- post.mode
proposal <- list(var=post.cov,scale=3)
set.seed(4234)
fit1.2 <- rwmetrop(logpost,proposal,start,T,y)
fit1.2$accept

# density plot
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2))
plot(density(fit1.2$par[5001:T,1]),main="",xlab=expression(theta[A]))
plot(density(fit1.2$par[5001:T,2]),main="",xlab=expression(theta[B]))

## (g)
# starting at the 1st mode
post.mode <- out1$par
post.cov <- -solve(out1$hessian)
post.sd <- sqrt(diag(post.cov))
set.seed(4234)
fit2.1 <- gibbs(logpost,start,T,2*post.sd,y)
fit2.1$accept

# density plots
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2))
plot(density(fit2.1$par[5001:T,1]),main="",xlab=expression(theta[A]))
plot(density(fit2.1$par[5001:T,2]),main="",xlab=expression(theta[B]))

# starting at the 2nd mode
post.mode <- out2$par
post.cov <- -solve(out2$hessian)
post.sd <- sqrt(diag(post.cov))
set.seed(4234)
fit2.2 <- gibbs(logpost,start,T,2*post.sd,y)
fit2.2$accept

# density plots
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2))
plot(density(fit2.2$par[5001:T,1]),main="",xlab=expression(theta[A]))
plot(density(fit2.2$par[5001:T,2]),main="",xlab=expression(theta[B]))


## (h) trace plots and autocorrelation plots
mcmcobj1.1 <- mcmc(fit1.1$par[5001:T,])
mcmcobj1.2 <- mcmc(fit1.2$par[5001:T,])
mcmcobj2.1 <- mcmc(fit2.1$par[5001:T,])
mcmcobj2.2 <- mcmc(fit2.2$par[5001:T,])
colnames(mcmcobj1.1) <- c("thetaA","thetaB")
colnames(mcmcobj1.2) <- c("thetaA","thetaB")
colnames(mcmcobj2.1) <- c("thetaA","thetaB")
colnames(mcmcobj2.2) <- c("thetaA","thetaB")

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(2,4))
traceplot(mcmcobj1.1)
traceplot(mcmcobj1.2)
autocorr.plot(mcmcobj1.1,auto.layout=FALSE)
autocorr.plot(mcmcobj1.2,auto.layout=FALSE)	

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(2,4))
traceplot(mcmcobj2.1)
traceplot(mcmcobj2.2)
autocorr.plot(mcmcobj2.1,auto.layout=FALSE)
autocorr.plot(mcmcobj2.2,auto.layout=FALSE)	
		
