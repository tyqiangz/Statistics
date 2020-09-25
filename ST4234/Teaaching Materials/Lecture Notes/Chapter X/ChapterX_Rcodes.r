###### Chapter X ######
###### The codes in this file include: ######
### 1. Bayesian linear regression
### 2. Bayesian model comparison 
### 3. Variational Bayes

###### 1. Bayesian linear regression
##### Example : Bird extinction data
require(invgamma)
require(coda)
require(mvtnorm)

require(LearnBayes)
data(birdextinct)
attach(birdextinct)
logtime=log(time)

# scatter plot of logtime vs nesting
par(mfrow=c(2,2))
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(nesting,logtime,col="blue")
out = (logtime > 3)
text(nesting[out], logtime[out], label=species[out], pos = 2)
plot(jitter(size),logtime,xaxp=c(0,1,1),col="blue")
plot(jitter(status),logtime,xaxp=c(0,1,1),col="blue")

# least square fit
fit1 <- lm(logtime~nesting+size+status,data=birdextinct,x=TRUE,y=TRUE)
summary(fit1)


## Bayesian inference
# find the posterior distribution and the predictive distribution by the following steps:
# 1. sampling from p(sigma2|y,X);
# 2. sampling from p(beta|sigma2,y,X);
# 3. predicting y

T <- 10^4
nu0 <- 1
sigma02 <- 0.01
g <- 100	# Zellner's g-prior

y <- logtime
X <- cbind(1,nesting,size,status)
n <- length(y)
p <- ncol(X)

(beta.ols <- solve(t(X)%*%X, t(X)%*%y))
beta.tilde <- g/(g+1) * beta.ols
S0 <- t(y)%*%y - t(beta.tilde)%*%((g+1)/g*t(X)%*%X)%*%beta.tilde

set.seed(4234)
sigma2.samples <- rinvgamma(T, shape=(nu0+n)/2, rate=(nu0*sigma02+S0)/2)
beta.samples <- matrix(0, nrow=T, ncol=p)
for(i in 1:T) {
	beta.samples[i,] <- rmvnorm(1, mean=beta.tilde, 
	                sigma=sigma2.samples[i]*g/(g+1)*solve(t(X)%*%X))
}

## posterior summary of beta and sigma2
apply(beta.samples,2,mean)
apply(beta.samples,2,sd)
apply(beta.samples,2,quantile,prob=c(.05,.5,.95))

mean(sigma2.samples)
sd(sigma2.samples)
quantile(sigma2.samples, prob=c(.05,.5,.95))

## predictive distribution for y at the observed locations 
error.pred <- matrix(0, nrow=T, ncol=n)
for (i in 1:T) {
	error.pred[i,] <- rnorm(n,mean=0,sd=sqrt(sigma2.samples[i]))
}
# y.pred are the draws from the predictive distribution
y.pred <- beta.samples%*%t(X) + error.pred
y.pred.CI <- matrix(0,nrow=n,ncol=2)
for (j in 1:n) {
	y.pred.CI[j,] <- quantile(y.pred[,j],prob=c(0.05,0.95))
}
     
# use plot of predictive CIs to find out the outliers
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
index <- 1:n
plot(index, y, pch=16, col="red",ylim=c(min(y.pred.CI)-0.5,max(y.pred.CI)+0.5))
for (i in 1:n){
	lines(rep(i,2),y.pred.CI[i,],lwd=2,lty=3,col="red")
	points(i,y.pred.CI[i,1],pch=3,col="blue")
	points(i,y.pred.CI[i,2],pch=3,col="blue")
}
outlier <- (y<y.pred.CI[,1] | y>y.pred.CI[,2])
text(index[outlier], y[outlier], label=species[outlier], pos = 4)


###### 2. Model Comparison
###### Example: Bayes factor 
###### Geometric model and Poisson model for count data
y <- c(0,0,0,0,0,0,1,1,1,1,1,2,2,3,4,4,4,5,6,7)
ysum <- sum(y)
n <- length(y)
m1 <- gamma(ysum+1)*gamma(n+1)/gamma(ysum+n+2)
m2 <- gamma(ysum+2)/(n+1)^(ysum+2)/prod(factorial(y))
m1
m2




###### 3. Variational Bayes
###### Example: Normal model
require(invgamma)

### Midge wing length (Hoff, 2009)
# input data
y <- c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)
n <- length(y)

# prior: mu|sigma2 ~ N(mu0, sigma2/n0), sigma2 ~ InvGamma(nu0/2, nu0*sigma02/2)
ybar <- mean(y)
s2 <- var(y)
mu0 <- 1.9
sigma02 <- 0.01
n0 <- 1
nu0 <- 1

mu1 <- (n*ybar+n0*mu0)/(n+n0)
S1 <- (n-1)*s2+nu0*sigma02+n*n0*(ybar-mu0)^2/(n+n0)

# VB distribution: q(mu)=N(u,v), q(sigma2)=InvGamma(a,b)
# parameters to be updated in VB distribution: v and b

# ELBO function
elbo.normal <- function(u,v,a,b,data) {
	n <- data$n
	n0 <- data$n0
	nu0 <- data$nu0
	sigma02 <- data$sigma02
	mu1 <- data$mu1
	S1 <- data$S1
	-(n+1)/2*log(2*pi)+nu0/2*log((nu0*sigma02)/2)-log(gamma(nu0/2))+1/2*log(n0)+
	((nu0+n+1)/2-a)*(digamma(a)-log(b))-
	1/2*(a/b)*((n+n0)*(u^2+v-2*u*mu1+mu1^2)+S1)+
	1/2+log(2*pi*v)/2-a*log(b)+log(gamma(a))+a
}

wingdata <- list(n=n,n0=n0,nu0=nu0,sigma02=sigma02,mu1=mu1,S1=S1)

# u and a are fixed
u <- mu1
a <- (nu0+n+1)/2

# initialize v and b
v <- 1
b <- 1
iter <- 0
elbo.diff <- 1
elbo.old <- elbo.normal(u,v,a,b,data=wingdata)

# find the mean field VB distribution using the CAVI algorithm
while(elbo.diff > 1e-6) {
	iter <- iter + 1
	v <- b/a/(n+n0)
	b <- ((n+n0)*(u^2+v-2*u*mu1+mu1^2)+S1)/2
	elbo.new <- elbo.normal(u,v,a,b,data=wingdata)
	elbo.diff <- abs(elbo.new - elbo.old)
	elbo.old <- elbo.new
	cat("iter=",iter,"; elbo=",elbo.new,"\n")
} 

u
v 
a
b
 
## contour plot
require(invgamma)
require(grDevices)

# true posterior function
joint.posterior <- function(para) {
	mu = para[1]
	sigma2 = para[2]
	dnorm(mu, mean=mu1, sd=sqrt(sigma2/(n+n0))) * 
	dinvgamma(sigma2, shape=(nu0+n)/2, rate=S1/2)
}

# mean field VB posterior function
mfvb.posterior <- function(para,u,v,a,b) {
	mu = para[1]
	sigma2 = para[2]
	dnorm(mu, mean=u, sd=sqrt(v)) * 
	dinvgamma(sigma2, shape=a, rate=b)
}

mu.grid <- seq(from=1.65, to=2, by=0.001)
sigma2.grid <- seq(from=0.001, to=0.04, by=0.001)
ts.grid <- expand.grid(mu.grid, sigma2.grid)
dens.true <- apply(ts.grid, 1, joint.posterior) 	# calculate true joint density
dens.vb <- apply(ts.grid, 1, mfvb.posterior, u=u, v=v, a=a, b=b) # calculate mean field VB density


par(mar=c(3.5,3.5,1.5,1))
par(mgp=c(2.1,0.8,0))
contour(mu.grid, sigma2.grid, matrix(dens.true,nrow=length(mu.grid),ncol=length(sigma2.grid)),
     col="blue",lwd=2, labcex=0.9, levels=seq(50,600,by=50),
	 xlab=expression(mu), ylab=expression(sigma^2), 
	 main="Contour plot of posteriors")	
contour(mu.grid, sigma2.grid, matrix(dens.vb,nrow=length(mu.grid),ncol=length(sigma2.grid)),
     col="red",lwd=2,lty=2,labcex=0.9, levels=seq(50,600,by=50), add=TRUE)	 
legend(1.9, 0.041, legend=c("True posterior", "VB posterior"),
       col=c("blue","red"), lwd=5,lty=c(1,2), bty="n")

	 
# calculate log p(y), the log marginal probability
(logpy <- log(gamma((nu0+n)/2))-log(gamma(nu0/2))+nu0/2*log(nu0*sigma02)-
         (nu0+n)/2*log(S1)-n/2*log(pi)+1/2*log(n0)-1/2*log(n+n0))

# verify that the calculation of log p(y) is correct
# set.seed(4234)
# s21 <- rinvgamma(10000,shape=nu0/2,rate=nu0*sigma02/2)
# theta1 <- rnorm(10000,mean=mu0,sd=sqrt(s21/n0))
# py <- vector()
# for(i in 1:10000){
	# py[i] <- prod(dnorm(y,mean=theta1[i],sd=sqrt(s21[i])))
# }
# log(mean(py))



