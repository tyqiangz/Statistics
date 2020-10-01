### Tutorial 3 codes

## Problem 3
# Poisson model with a Gamma prior
y <- c(0,0,0,0,0,0,1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,5)
n <- length(y)	   

# (a) Gamma(2,1) prior	   
a0 <- 2		# prior parameters
b0 <- 1

a <- a0 + sum(y)	# posterior parameters
b <- b0 + n

theta <- seq(from=0, to=5, by=0.001)
plot(theta, dgamma(theta,shape=a0,rate=b0), type="l", lwd=2.5,
     xlab=expression(theta), ylab="Density", ylim=c(0,2.2))
points(theta, dgamma(theta,shape=a,rate=b), type="l", lwd=2.5, col="red")
legend(3,2,legend = c(expression(paste("p(", theta, ")")),
       expression(paste("p(", theta, "|y)"))), lty=1, 
       lwd=2.5, cex=1, col=c("black","red"), bty="n")

# posterior mean
a/b
# posterior standard deviation
sqrt(a/b^2)
# posterior mode
(a-1)/b

# (b) Gamma(1/2,0) prior
a0 <- 1/2		# prior parameters
b0 <- 0

a <- a0 + sum(y)	# posterior parameters
b <- b0 + n

theta <- seq(from=0, to=5, by=0.001)
plot(theta, theta^(-1/2), type="l", lwd=2.5,
     xlab=expression(theta), ylab="Density", ylim=c(0,2.2))
points(theta, dgamma(theta,shape=a,rate=b), type="l", lwd=2.5, col="red")
legend(3,2,legend = c(expression(paste(theta^"-1/2")),
       expression(paste("p(", theta, "|y)"))), lty=1, 
       lwd=2.5, cex=1, col=c("black","red"), bty="n")

# posterior mean
a/b
# posterior standard deviation
sqrt(a/b^2)
# posterior mode
(a-1)/b

	   