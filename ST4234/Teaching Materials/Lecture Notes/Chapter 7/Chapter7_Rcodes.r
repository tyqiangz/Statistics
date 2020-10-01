###### Chapter 7 MCMC: Metropolis-Hastings Algorithm ######

### simulate a discrete Markov chain (Page 8)

# transition matrix
P <- rbind(c(0.5,0.5,0,0,0,0),
           c(0.25,0.5,0.25,0,0,0),
           c(0,0.25,0.5,0.25,0,0),
           c(0,0,0.25,0.5,0.25,0),
           c(0,0,0,0.25,0.5,0.25),
           c(0,0,0,0,0.5,0.5))

N <- 50000         # total no. of steps
l <- rep(0,N)      # path taken 
l[1] <- 3          # starting location = 3

set.seed(4234)
for (i in 2:N){
	l[i] <- sample(1:6, size=1, prob = P[l[i-1],])
}

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(l[1:200],type="l",xlab="step",ylab="state",col="blue",
     cex.axis=1.2,cex.lab=1.4)

table(l[1:500])/500
table(l[1:5000])/5000
table(l[1:50000])/50000

w <- matrix(c(0.1,0.2,0.2,0.2,0.2,0.1),nrow=1)
w%*%P



### Example: Learning about normal population from grouped data (Page 40)
data <- cbind( c(-Inf,66,68,70,72,74),
               c(66,68,70,72,74,Inf),
               c(14,30,49,70,33,15) )

# compute log posterior
dpnorm <- function(a,b,mu,sigma){
	pnorm(b,mean=mu,sd=sigma) - pnorm(a,mean=mu,sd=sigma)
}

logpost <- function(theta, data){
	G <- nrow(data)
	L <- 0
	mu <- theta[1]
	sigma <- exp(theta[2])
	for (g in 1:G){
		L <- L + data[g,3]*log(dpnorm(data[g,1],data[g,2],mu,sigma))
	}
	return(L)
}

# contour plot
N <- 100
theta1 <- seq(from=66, to=74, length.out=N)
theta2 <- seq(from=0, to=2.5, length.out=N)
z <- matrix(0,N,N)
for (i in 1:N){
	for (j in 1:N){
		theta <- c(theta1[i], theta2[j])
		z[i,j] <- logpost(theta,data)    
}}

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
contour(x=theta1,y=theta2,z,col="blue",nlevels=30,
        xlab=expression(mu),ylab=expression(lambda),
        cex.axis=1.1,cex.lab=1.3,labcex=1.1)

# find normal approximation
start <- c(70,1)
out <- optim(par=start,fn=logpost,hessian=TRUE,
             control=list(fnscale=-1),data=data)
(post.mode <- out$par)
(post.cov <- -solve(out$hessian))


# Random walk Metropolis-Hastings
require(LearnBayes)
proposal <- list(var=post.cov,scale=2)
set.seed(4234)
fit1 <- rwmetrop(logpost,proposal,start,10000,data)
fit1$accept
# you can also start at post.mode directly; it converges faster

(post.means <- apply(fit1$par[5001:10000,],2,mean))
(post.sds <- apply(fit1$par[5001:10000,],2,sd))

post.mode
sqrt(diag(post.cov))


# contour plot of the posterior (together with the MH draws)
N <- 100
theta1 <- seq(from=69, to=71.4, length.out=N)
theta2 <- seq(from=0.65, to=1.3, length.out=N)
z <- matrix(0,N,N)
for (i in 1:N){
	for (j in 1:N){
		theta <- c(theta1[i], theta2[j])
		z[i,j] <- logpost(theta,data)    
}}

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
contour(x=theta1,y=theta2,z,col="blue",nlevels=20,
        xlab=expression(mu),ylab=expression(lambda),
        cex.axis=1.1,cex.lab=1.3,labcex=1.1)
points(fit1$par[5001:10000,])


# output analysis
start <- c(65,1)
proposal <- list(var=post.cov,scale=0.2)
# "Algorithm 2"
set.seed(4234)
fit2 <- rwmetrop(logpost,proposal,start,10000,data)
fit2$accept

# create mcmc objects (if input is a matrix, 
# each column should represent one variable)
require(coda)
colnames(fit1$par) <- c("mu","lambda")
colnames(fit2$par) <- c("mu","lambda")
mcmcobj1 <- mcmc(fit1$par)
mcmcobj2 <- mcmc(fit2$par)

# traceplots
par(mar=c(3.5,2,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2),oma=c(0,0,1,0))
traceplot(mcmcobj1)
title("Algorithm 1", outer=TRUE)

par(mar=c(3.5,2,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2),oma=c(0,0,1,0))
traceplot(mcmcobj2)
title("Algorithm 2", outer=TRUE)

# remove burn-in period
mcmcobj1 <- mcmc(fit1$par[2001:10000,])
mcmcobj2 <- mcmc(fit2$par[2001:10000,])

# traceplots
par(mar=c(3.5,2,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2),oma=c(0,0,1,0))
traceplot(mcmcobj1,cex.lab=1.2)
title("Algorithm 1", outer=TRUE, cex.main=1.5)

par(mar=c(3.5,2,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2),oma=c(0,0,1,0))
traceplot(mcmcobj2,cex.lab=1.2)
title("Algorithm 2", outer=TRUE, cex.main=1.5)

# autocorrelation plots
par(mar=c(3.5,2,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2),oma=c(0,0,1,0))
autocorr.plot(mcmcobj1)
title("Algorithm 1", outer=TRUE, cex.main=1.5)

par(mar=c(3.5,2,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2),oma=c(0,0,1,0))
autocorr.plot(mcmcobj2)
title("Algorithm 2", outer=TRUE, cex.main=1.5)


# summary statistics
summary(mcmcobj1)
summary(mcmcobj2)



### A direct R implementation of MH algorithm for the same model
# the same logpost() function as before will be used
require(mvtnorm)
T <- 10^4
theta.trace <- matrix(0,nrow=T,ncol=2)	# record the MCMC draws of theta=(mu,lambda)
colnames(theta.trace) <- c("mu","lambda")
theta.old <- c(70,1)	# initial value
logpost.old <- logpost(theta.old, data) # current value of log posterior
scale.V <- 4*post.cov	# proposal covariance matrix; 
# note that scale=2 means you need to multiply post.cov by 2^2=4
accept <- numeric(T)	# record acceptance=1, rejection=0, for each iteration

set.seed(4234)
#system.time(
for (i in 1:T){
	theta.new <- as.vector(rmvnorm(1, mean=theta.old, sigma=scale.V))
	logpost.new <- logpost(theta.new, data)
	alpha <- min(1, exp(logpost.new - logpost.old))
	u = runif(1)
	if (u < alpha) {
		accept[i] = 1
		theta.old = theta.new
		logpost.old = logpost.new	# update log likelihood only if theta changes
	} else {
		accept[i] = 0
	}
	theta.trace[i,] = theta.old
}
#)
# compare with system.time(fit1 <- rwmetrop(logpost,proposal,start,10000,data))

# trace plot of mu and lambda
mcmcobj3 <- mcmc(theta.trace)

# traceplots
par(mar=c(3.5,2,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2),oma=c(0,0,1,0))
traceplot(mcmcobj3)
title("Direct Coding MH", outer=TRUE)

# posterior mean and sd of mu and lambda
(post.means <- apply(theta.trace[5001:T,],2,mean))
(post.sds <- apply(theta.trace[5001:T,],2,sd))

# acceptance rate
mean(accept)



  