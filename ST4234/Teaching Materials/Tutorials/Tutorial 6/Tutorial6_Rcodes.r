### Tutorial 6 codes

## Problem 2

# define prior hyperparameters
mu <- 0
sigma2 <- 0.25^2
n <- 5
y <- 5	# 5 heads

# log posterior function
logpost <- function(theta,y,n,mu,sigma2){
	theta*y - n*log(1+exp(theta)) - 0.5*(theta-mu)^2/sigma2
}

# plot the log posterior
theta.grid <- seq(from=-10,to=10,by=0.1)
plot(theta.grid, logpost(theta.grid,y,n,mu,sigma2),type="l",lwd=2,
     ylab=expression(paste("h(",theta,")")),xlab=expression(theta),
	 main="log posterior")

	 
# normal approximation: method 1, using optimize()
(out <- optimize(logpost, interval=c(-10,10), maximum=TRUE, y=y,n=n,mu=mu,sigma2=sigma2))
(post.mode <- out$maximum)	 
(post.var <- 1/(n*exp(post.mode)/(1+exp(post.mode))^2 + 1/sigma2))	# variance using analytical Hessian 
	 
	 
# normal approximation: method 2, using optim() to find the mode and Hessian
(out <- optim(par=0, fn=logpost, hessian=TRUE, control=list(fnscale=-1),
        method="Brent", lower=-10, upper=10,y=y,n=n,mu=mu,sigma2=sigma2))
(post.mode <- out$par)
(post.var <- -1/out$hessian)

# find P(theta>0|y=5) using normal approximation 
1 - pnorm(0, mean=post.mode, sd=sqrt(post.var))
	 
	 
	 
## Problem 3
# (b) log posterior function
s <- 8
n <- 15
t <- 15962989
t1 <- 237217
logpost <- function(theta, s, n, t, t1){
	theta1 <- theta[1]
	theta2 <- theta[2]
	return((1-s)*theta1 + theta2 - (t-n*t1)*exp(-theta1) - n*exp(theta2-theta1))
}

# (c) 
# some trial-and-error / exploration
(out <- optim(par=c(0,0), fn=logpost, hessian=TRUE,
              control=list(fnscale=-1), s=s, n=n, t=t, t1=t1))
# this ends up in a local extremum value, not the global maximum
			  
# contour plot: several attempts
# theta1.grid <- seq(from=-100, to=100, by=1)
# theta2.grid <- seq(from=-100, to=100, by=1)	# range too large, change the range
theta1.grid <- seq(from=13, to=20, by=.1)
theta2.grid <- seq(from=-10, to=18, by=.1)
theta.grid <- expand.grid(theta1.grid, theta2.grid)			  
grid.logpost <- apply(theta.grid, 1, logpost, s=s, n=n, t=t, t1=t1) 			  
contour(theta1.grid, theta2.grid, matrix(grid.logpost,nrow=length(theta1.grid),ncol=length(theta2.grid)),
        col="blue", nlevels=800, lwd=2, labcex=0.9, 
		xlab=expression(theta[1]), ylab=expression(theta[2]))				  

# change the starting point to (14.5,11)
(out <- optim(par=c(14.5,11), fn=logpost, hessian=TRUE,
              control=list(fnscale=-1), s=s, n=n, t=t, t1=t1))		
(post.mode <- out$par)
(post.var <- -solve(out$hessian))		
		