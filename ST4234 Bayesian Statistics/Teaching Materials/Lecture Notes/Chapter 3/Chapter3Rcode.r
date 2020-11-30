###### Chapter 3: Normal model	######


### Example: Midge wing length (Hoff, 2009)
# input data
y <- c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)
n <- length(y)

## prior specification, when theta is unknown and sigma^2 is known
# prior of theta: N(mu0, tau02)
ybar <- mean(y)
sigma2 <- var(y)
mu0 <- 1.9
tau02 <- 0.95^2

# posterior parameters
# posterior of theta: N(mun, taun2)
(taun2 <- 1 / (n/sigma2 + 1/tau02))
(mun <- taun2 * (n*ybar/sigma2 + mu0/tau02))

# 95% Bayesian CI
qnorm(c(0.025,0.975), mean=mun, sd=sqrt(taun2))


## prior specification, when both theta and sigma^2 are unknown
# prior: theta|sigma2 ~ N(mu0, sigma2/n0), sigma2 ~ InvGamma(nu0/2, nu0*sigma02/2)
s2 <- var(y)
mu0 <- 1.9
sigma02 <- 0.01
n0 <- 1
nu0 <- 1

# posterior parameters
# posterior: theta|sigma2,y ~ N(mu1, sigma2/n1), sigma2|y ~ InvGamma(nu1/2, nu1*sigma12/2)
(n1 <- n + n0)
(mu1 <- (n*ybar+n0*mu0)/n1)
(nu1 <- n + nu0)
(sigma12 <- (nu0*sigma02 + (n-1)*s2 + n*n0*(ybar-mu0)^2/n1)/nu1)

# joint posterior of (theta,sigma2)
require(invgamma)
require(grDevices)
joint.posterior <- function(para) {
	theta = para[1]
	sigma2 = para[2]
	dnorm(theta, mean=mu1, sd=sqrt(sigma2/n1)) * 
	dinvgamma(sigma2, shape=nu1/2, rate=nu1*sigma12/2)
}

theta.grid <- seq(from=1.6, to=2, by=0.001)
sigma2.grid <- seq(from=0.001, to=0.04, by=0.001)
ts.grid <- expand.grid(theta.grid, sigma2.grid)
dens <- apply(ts.grid, 1, joint.posterior) 	# calculate joint density

# Plots
par(mfrow=c(1,2))
par(mar=c(3.5,3.5,1.5,1))
par(mgp=c(2.1,0.8,0))
# Contour plot
contour(theta.grid, sigma2.grid, matrix(dens,nrow=length(theta.grid),ncol=length(sigma2.grid)),
     col="blue",lwd=2, labcex=0.9, xlab=expression(theta), ylab=expression(sigma^2),
	 main="Contour plot of joint posterior")	 
# Image plot
image(theta.grid, sigma2.grid, matrix(dens,nrow=length(theta.grid),ncol=length(sigma2.grid)),
	  xlab=expression(theta), ylab=expression(sigma^2), main="Image plot of joint posterior")



