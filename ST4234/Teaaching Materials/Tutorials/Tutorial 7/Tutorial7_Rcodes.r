###### Tutorial 7 codes

### Problem 1
## (b)
# log posterior function
logpost <- function(eta) {
	125*log(2+(exp(eta)/(1+exp(eta))))-74*log(1+exp(eta))+35*eta	
}

# plot the log posterior
eta.grid <- seq(from=-10, to=10, by=0.1)
plot(eta.grid, logpost(eta.grid),type="l",lwd=2,
     ylab=expression(paste("logpost(",eta,")")),xlab=expression(eta),
	 main=expression(paste("log posterior of ",eta)))

# find the posterior mode using optim() with method="Brent"
(out <- optim(par=0, fn=logpost, hessian=TRUE, control=list(fnscale=-1), 
              method="Brent", lower=-10, upper=10))
(post.mode <- out$par)
(post.var <- -1/out$hessian)

# normal 95% CI for eta
(CI <- qnorm(c(0.025,0.975),mean=post.mode,sd=sqrt(post.var)))
# transformed 95% CI for theta
exp(CI)/(1+exp(CI))			  

	
## (c) Rejection sampling

# difference function = log(f) - log(g)
require(mvtnorm)
df <- 4
diff <- function(eta, post.mode, post.var, df){
    logpost(eta) - dmvt(matrix(eta,ncol=1), delta=post.mode, 
	sigma=matrix(post.var,1,1), df=df)
}

# maximize the difference function
(out <- optim(par=0, fn=diff, control=list(fnscale=-1),
              method="Brent", lower=-10, upper=10,
              post.mode=post.mode, post.var=post.var, df=df))
logM <- out$value	# M' value

# rejection sampling function
rej_sampling <- function(S, post.mode, post.var, df, logM){
	eta.samples <- rep(0,S)
	s <- 0
	T <- 0
	while (s < S){
		T <- T+1
		eta <- rmvt(1, delta=post.mode, sigma=matrix(post.var,1,1), df=df)
		U <- runif(1)
		if (log(U) < diff(eta, post.mode, post.var, df) - logM){
			s <- s+1
			eta.samples[s] <- eta
		}
	}
	accept.rate <- S/T
	list(accept.rate=accept.rate, eta.samples=eta.samples)
}			  

set.seed(4234)
S <- 10^4
output <- rej_sampling(S, post.mode, post.var, df, logM)
output$accept.rate

# CI for eta
quantile(output$eta.samples, c(0.025,0.975))

# CI for theta
theta.samples <- 1/(1+exp(-output$eta.samples))
quantile(theta.samples, c(0.025,0.975))

	

### Problem 2
## (a) 
# first, define the log posterior function
logpost <- function(theta, Y) {
	Y[1]*log(2+theta) + (Y[2]+Y[3])*log(1-theta) + Y[4]*log(theta)
}

# Y=(125,18,20,34)
Y <- c(125,18,20,34)
# optimize logpost()
# find the posterior mode and variance of normal approximation
# formulas can be found in Chapter 5 notes
(out <- optimize(logpost, interval=c(0,1), Y=Y, maximum=TRUE))
post.mode <- out$maximum
post.var <- 1/(Y[1]/(2+post.mode)^2+(Y[2]+Y[3])/(1-post.mode)^2+Y[4]/post.mode^2)
# importance sampling: draw samples from normal
set.seed(4234)
S <- 10^4
theta.samples <- rnorm(S, mean=post.mode, sd=sqrt(post.var))
# raw log weights
logw <- logpost(theta.samples,Y=Y) - 
        dnorm(theta.samples, mean=post.mode, sd=sqrt(post.var), log=TRUE)
w <- exp(logw - max(logw))	# raw weights
W <- w/sum(w)		# normalized weights
hist(W, breaks=30)		# histogram of normalized weights
(est <- sum(W*theta.samples))		# IS estimate of E(theta|y)
(se <- sqrt(sum((W*(theta.samples -est))^2)))	# standard error
		
# Y=(14,0,1,5)
Y <- c(14,0,1,5)
# optimize logpost()
# find the posterior mode and variance of normal approximation
# formulas can be found in Chapter 5 notes
(out <- optimize(logpost, interval=c(0,1), Y=Y, maximum=TRUE))
post.mode <- out$maximum
post.var <- 1/(Y[1]/(2+post.mode)^2+(Y[2]+Y[3])/(1-post.mode)^2+Y[4]/post.mode^2)
# importance sampling: draw samples from normal
set.seed(4234)
S <- 10^4
theta.samples <- rnorm(S, mean=post.mode, sd=sqrt(post.var))
# raw log weights
logw <- logpost(theta.samples,Y=Y) - 
        dnorm(theta.samples, mean=post.mode, sd=sqrt(post.var), log=TRUE)
# logw contains NaNs for those theta's falling outside [0,1]
# remove all NaNs from logw
# find the indices where logw is not NaN, i.e. theta is inside [0,1]
index <- which(is.nan(logw)==FALSE)
(length(index))		# number of theta samples inside [0,1]
w <- exp(logw[index] - max(logw[index]))	# raw weights
W <- w/sum(w)		# normalized weights
hist(W, breaks=30)		# histogram of normalized weights
(est <- sum(W*theta.samples[index]))		# IS estimate of E(theta|y)
(se <- sqrt(sum((W*(theta.samples[index] -est))^2)))	# standard error	


## (b)
# Y=(125,18,20,34)
Y <- c(125,18,20,34)
(out <- optimize(logpost, interval=c(0,1), Y=Y, maximum=TRUE))
post.mode <- out$maximum
post.var <- 1/(Y[1]/(2+post.mode)^2+(Y[2]+Y[3])/(1-post.mode)^2+Y[4]/post.mode^2)
# find beta parameters: define a 2d loss function
beta.para <- function(para, post.mode, post.var) {
	a <- para[1]; b <- para[2]
	ab.mode = (a-1)/(a+b-2)
	ab.var = a*b/(a+b)^2/(a+b+1)
	return((ab.mode-post.mode)^2+(ab.var-post.var)^2)
}
(out.ab <- optim(par=c(5,5), fn=beta.para, 
                post.mode=post.mode, post.var=post.var))
a <- out.ab$par[1]
b <- out.ab$par[2]
# importance sampling: draw samples from Beta(a,b)
set.seed(4234)
S <- 10^4
theta.samples <- rbeta(S, a, b)
# raw log weights
logw <- logpost(theta.samples,Y=Y) - 
        dbeta(theta.samples, a, b, log=TRUE)
w <- exp(logw - max(logw))	# raw weights
W <- w/sum(w)		# normalized weights
hist(W, breaks=30)		# histogram of normalized weights
(est <- sum(W*theta.samples))		# IS estimate of E(theta|y)
(se <- sqrt(sum((W*(theta.samples -est))^2)))	# standard error


# Y=(14,0,1,5)
Y <- c(14,0,1,5)
(out <- optimize(logpost, interval=c(0,1), Y=Y, maximum=TRUE))
post.mode <- out$maximum
post.var <- 1/(Y[1]/(2+post.mode)^2+(Y[2]+Y[3])/(1-post.mode)^2+Y[4]/post.mode^2)
# find beta parameters: define a 2d loss function
beta.para <- function(para, post.mode, post.var) {
	a <- para[1]; b <- para[2]
	ab.mode = (a-1)/(a+b-2)
	ab.var = a*b/(a+b)^2/(a+b+1)
	return((ab.mode-post.mode)^2+(ab.var-post.var)^2)
}
(out.ab <- optim(par=c(5,5), fn=beta.para, 
                 post.mode=post.mode, post.var=post.var))
a <- out.ab$par[1]
b <- out.ab$par[2]
# importance sampling: draw samples from Beta(a,b)
set.seed(4234)
S <- 10^4
theta.samples <- rbeta(S, a, b)
# raw log weights
logw <- logpost(theta.samples,Y=Y) - 
        dbeta(theta.samples, a, b, log=TRUE)
w <- exp(logw - max(logw))	# raw weights
W <- w/sum(w)		# normalized weights
hist(W, breaks=30)		# histogram of normalized weights
(est <- sum(W*theta.samples))		# IS estimate of E(theta|y)
(se <- sqrt(sum((W*(theta.samples-est))^2)))	# standard error



	
### Problem 3
mu <- 0
sigma2 <- 0.25^2
n <- 5
y <- 5
logpost <- function(theta, y, n, mu, sigma2){
	theta*y - n*log(1+exp(theta)) - 0.5*(theta-mu)^2/sigma2
}
diff <- function(theta,y,n,mu,sigma2){
	logpost(theta,y,n,mu,sigma2) - 
	dnorm(theta,mean=mu,sd=sqrt(sigma2),log=TRUE)
}

# find M' by maximizing the diffence function
(out <- optim(par=10, fn=diff,control=list(fnscale=-1),method="Brent",
              lower=0,upper=100,y=y,n=n,mu=mu,sigma2=sigma2))
logM <- out$value

# plot of the difference function
theta.grid <- seq(from=-20,to=50,by=0.1)
plot(theta.grid, diff(theta.grid,y,n,mu,sigma2),type="l",
	ylab=expression(paste("logpost(",theta,") - log(p(",theta,"))")),
	xlab=expression(theta))
# Alternatively, by analytical calculation
(logM <- 0.5*log(2*pi)+0.5*log(sigma2))

# rejection sampling function
rej_sampling <- function(S, y, n, mu, sigma2, logM){
	theta.samples <- rep(0,S)
	s <- 0
	T <- 0
	while (s < S){
		T <- T+1
		theta <- rnorm(1, mean=mu, sd=sqrt(sigma2))
		U <- runif(1)
		if (log(U) < diff(theta, y, n, mu, sigma2) - logM){
			s <- s+1
			theta.samples[s] <- theta
		}
	}
	accept.rate <- S/T
	list(accept.rate=accept.rate, theta.samples=theta.samples)
}				  
			  
set.seed(4234)
S <- 10^4
output <- rej_sampling(S, y, n, mu, sigma2, logM)
output$accept.rate	

# estimate P(theta>0|y)
mean(output$theta.samples>0)




		  