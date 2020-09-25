### Tutorial 2 codes

## Problem 1
# (a) 95% quantile-based CI
qbeta(c(0.025,0.975), 17, 36)

# (a) Plots of prior, likelihood, and posterior
theta <- seq(from=0,to=1,by=0.001)
y <- 15;  n <- 43
a0 <- 2;  b0 <- 8
a <- a0 + y;  b <- b0 + n - y
likelihood <- dbinom(y, n, theta)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2))
plot(theta, likelihood, type="l", lwd=2.5, col="blue",
     xlab=expression(theta), ylab="Density")
legend(0.45,0.12,legend=expression(paste("p(y|", theta, ")")),
       lty=1, lwd=2.5, cex=1, col="blue", bty="n")
plot(theta, dbeta(theta,a0,b0), type="l", lwd=2.5,
     xlab=expression(theta), ylab="Density", ylim=c(0,6.2))
points(theta, dbeta(theta,a,b), type="l", lwd=2.5, col="red")
legend(0.45,6,legend = c(expression(paste("p(", theta, ")")),
       expression(paste("p(", theta, "|y)"))), lty=1, 
       lwd=2.5, cex=1, col=c("black","red"), bty="n")

	   
	   
# (b) 95% quantile-based CI
qbeta(c(0.025,0.975), 23, 30)

# (b) Plots of prior, likelihood, and posterior
theta <- seq(from=0,to=1,by=0.001)
y <- 15;  n <- 43
a0 <- 8;  b0 <- 2
a <- a0 + y;  b <- b0 + n - y
likelihood <- dbinom(y, n, theta)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2))
plot(theta, likelihood, type="l", lwd=2.5, col="blue",
     xlab=expression(theta), ylab="Density")
legend(0.45,0.12,legend=expression(paste("p(y|", theta, ")")),
       lty=1, lwd=2.5, cex=1, col="blue", bty="n")
plot(theta, dbeta(theta,a0,b0), type="l", lwd=2.5,
     xlab=expression(theta), ylab="Density", ylim=c(0,6.2))
points(theta, dbeta(theta,a,b), type="l", lwd=2.5, col="red")
legend(0.45,6,legend = c(expression(paste("p(", theta, ")")),
       expression(paste("p(", theta, "|y)"))), lty=1, 
       lwd=2.5, cex=1, col=c("black","red"), bty="n")
	   
	   
# (c) plots of three priors
beta_mixture <- function(theta){
	(3/4)*dbeta(theta,2,8) + (1/4)*dbeta(theta,8,2)
}
theta <- seq(from=0, to=1, by=0.001)
plot(theta, beta_mixture(theta), type="l", lwd=2.5,
     ylim=c(0,4), xlab=expression(theta), ylab="Density",
	 main="Prior densities of theta")
points(theta, dbeta(theta,2,8), type="l", lwd=2.5, col="red")
points(theta, dbeta(theta,8,2), type="l", lwd=2.5, col="blue")
legend(0.15, 4.2, legend = c(expression(paste("p"[1], "(", theta, ")")),
       expression(paste("p"[2], "(", theta, ")")),
	   expression(paste("p(",theta,")=(3/4)p"[1],"(",theta,")",
	   "+(1/4)p"[2],"(",theta,")"))), 
	   lty=1, lwd=2.5, cex=1, bty="n", col=c("red","blue","black"))

# (d)-(iii) Plot of the posterior. Find the mode.
y <- 15
n <- 43
posterior_prop <- function(theta,y,n){
	beta_mixture(theta) * dbinom(y,n,theta)
}
plot(theta, posterior_prop(theta,y=y,n=n), type="l", lwd=2.5,
	xlab=expression(theta), ylab="Density",
	main=expression(paste("p(", theta, ") x p(y|", theta, ")"))) 	 
(posterior_mode <- optimize(posterior_prop, interval=c(0,1),
 maximum=TRUE, n=n, y=y))
	

## Problem 2
# (b) Plots of Galenshore family.
theta <- seq(from=0,to=5, by=0.001)
dprior <- function(theta,tau0,nu0){
2*tau0^(2*nu0)/gamma(nu0)*theta^(2*nu0-1)*exp(-tau0^2*theta^2)
}
plot(theta, dprior(theta,0.1,0.1), type="l",lwd=2.5,
xlab=expression(theta), ylab="Density",
main="Galenshore priors",ylim=c(0,5))
points(theta, dprior(theta,1,1), type="l",lwd=2.5,col="blue")
points(theta, dprior(theta,3,3), type="l",lwd=2.5,col="green")
points(theta, dprior(theta,5,1), type="l",lwd=2.5,col="red")
points(theta, dprior(theta,1,5), type="l",lwd=2.5,col="orange")
legend("topright",
legend=c(expression(paste(tau[0]," = ",nu[0]," = 0.1")),
expression(paste(tau[0]," = ",nu[0]," = 1")),
expression(paste(tau[0]," = ",nu[0]," = 3")),
expression(paste(tau[0]," = 5,",nu[0]," = 1")),
expression(paste(tau[0]," = 1,",nu[0]," = 5"))),
col=c("black","blue","green","red","orange"),
lty=1,lwd=2.5)
	
	
