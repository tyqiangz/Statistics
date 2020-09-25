###### Chapter 4 Monte Carlo Approximation	######

### Monte Carlo approximation of gamma distribution

set.seed(1)   
Svalues <- c(10,100,1000)
theta <- seq(from=0.5, to=2.5, by=0.01) 
par(mfcol=c(2,3))
par(mar=c(3.5,3.5,1,1))
par(mgp=c(1.9,0.7,0))
for (i in 1:3){
	S <- Svalues[i]
	samples <- rgamma(S,68,45)
	hist(samples,freq=FALSE,xlim=c(0.5,2.5),ylim=c(0,2.5),
         xlab=expression(theta),main=paste("S =",S))
	points(theta,dgamma(theta,68,45),type="l",col="red")  
	plot(density(samples),xlim=c(0.5,2.5),ylim=c(0,2.5),
         xlab=expression(theta),main=paste("S =",S))
	points(theta,dgamma(theta,68,45),type="l",col="red")  
}


# Expectation
set.seed(1)
Svalues <- c(10,100,1000)
m <- rep(0,3)
for (i in 1:3){
	S <- Svalues[i]
	samples <- rgamma(S,68,45)
	m[i] <- mean(samples)
}
m


# Probabilities
pgamma(1.75,68,45)
set.seed(1)
Svalues <- c(10,100,1000)
prob <- rep(0,3)
for (i in 1:3){
	S <- Svalues[i]
	samples <- rgamma(S,68,45)
	prob[i] <- sum(samples<1.75)/S
}
prob


# Confidence intervals
qgamma(c(0.025,0.975),68,45)
set.seed(1)
Svalues <- c(10,100,1000)
CI <- matrix(0,3,2)
for (i in 1:3){
	S <- Svalues[i]
	samples <- rgamma(S,68,45)
	CI[i,] <- quantile(samples, c(0.025,0.975))
}
CI


# Convergence (with increasing sample sizes, used for the plot on page 13)
set.seed(1)
m <- rep(0,1000)
prob <- rep(0,1000)
upperCI <- rep(0,1000)
samples <- rgamma(1000,68,45)
for (i in 1:1000){
	m[i] <- mean(samples[1:i])
	prob[i] <- sum(samples[1:i]<1.75)/i
	upperCI[i] <- quantile(samples[1:i], 0.975)
}


# plots of cumulative estimates on page 13
par(mfrow=c(1,3))
par(mar=c(3.5,3.5,1,1))
par(mgp=c(1.9,0.7,0))
plot(m,type="l",xlab="No. of Monte Carlo samples", ylab="cumulative mean",cex.lab=1.2,cex.axis=1.2)
abline(h=68/45,col="red",lwd=2)
plot(prob,type="l",xlab="No. of Monte Carlo samples", ylab="cumulative cdf at 1.75",cex.lab=1.2,cex.axis=1.2)
abline(h=pgamma(1.75,68,45),col="red",lwd=2)
plot(upperCI,type="l",xlab="No. of Monte Carlo samples", ylab="cumulative 97.5% quantile",cex.lab=1.2,cex.axis=1.2)  
abline(h=qgamma(0.975,68,45),col="red",lwd=2)


### Example: log-odds
set.seed(1)
S <- 10000
a0 <- 1; b0 <- 1
y <- 441; n <- 860
theta_samples_prior <- rbeta(S,a0,b0)
gamma_samples_prior <- log(theta_samples_prior/(1-theta_samples_prior))

theta_samples_post <- rbeta(S,a0+y,b0+n-y)
gamma_samples_post <- log(theta_samples_post/(1-theta_samples_post))


par(mfrow=c(1,2))
par(mar=c(3.5,3.5,1,1))
par(mgp=c(1.9,0.7,0))
plot(density(gamma_samples_prior), main="prior", xlab=expression(gamma),xlim=c(-4.5,4.5)) 
plot(density(gamma_samples_post), main="posterior", xlab=expression(gamma),xlim=c(-4.5,4.5))


### Example: Birth rates
# posterior probability of theta1>theta2
set.seed(1)
S <- 10^4
theta1.samples <- rgamma(S,219,112)
theta2.samples <- rgamma(S,68,45)
sum(theta1.samples>theta2.samples)/S

par(mar=c(3.5,3.5,1,1))
par(mgp=c(1.9,0.7,0))
gamma.samples <- theta1.samples/theta2.samples
plot(density(gamma.samples), xlab=expression(gamma),main="",cex.lab=1.3)


# predictive distribution

S <- 10^6

# Method 1: Use the exact predictive distribution, which is negative binomial
# note: the exact predictive distribution can be found on page 59 of Chapter 2
# note: this method is discussed on page 34 of Chapter 4 notes
set.seed(1)
r1 <- 219;  p1 <- 112/113
r2 <- 68;   p2 <- 45/46
predy1.samples <- rnbinom(S, r1, p1)
predy2.samples <- rnbinom(S, r2, p2)
sum(predy1.samples > predy2.samples)/S


# Method 2: first draw from the posterior (gamma distribution), then draw from Poisson
# note: this method is discussed on pages 32-33 of Chapter 4 notes
set.seed(1)
theta1.samples <- rgamma(S,219,112)
theta2.samples <- rgamma(S,68,45)
predy1.samples <- rpois(S,theta1.samples)
predy2.samples <- rpois(S,theta2.samples)
sum(predy1.samples > predy2.samples)/S


# summarize and plot the distribution of difference D = Y1 - Y2
D <- predy1.samples - predy2.samples
summary(D)
D <- as.factor(D)
summary(D)
prob <- summary(D)/S
sum(prob)

par(mar=c(3.5,3.5,1,1))
par(mgp=c(1.9,0.7,0))
plot(as.numeric(levels(D)), prob, type="h", 
    xlab="D",ylab="density",lwd=2,cex.lab=1.2,cex.axis=1.2)


# women without college degrees

nchildren <- 0:9
nwomen <- c(20,19,38,20,10,2,2,0,0,0)
n1 <- 111
par(mar=c(3.5,3.5,1,1))
par(mgp=c(1.9,0.7,0))
plot(nchildren-0.07,nwomen/n1,type="h",lwd=5,xlab="no. of children",ylab="probability",cex.lab=1.2,cex.axis=1.2)
points(nchildren+0.07,dnbinom(nchildren,219,112/113),type="h",col="red",lwd=5)
legend("topright", legend=c("empirical distribution", "predictive distribution"),
    col=c("black","red"),lty=1,lwd=5,cex=1.1)


set.seed(1)
S <- 10^4
n1 <- 111
theta.samples <- rgamma(S,219,112)
t.samples <- rep(0,S)
for (s in 1:S){
  Y.samples <- rpois(n1,theta.samples[s])
  t.samples[s] <- sum(Y.samples==2)/sum(Y.samples==1)
}

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2,0.7,0))
plot(density(t.samples),xlab="t",ylab="density",main="",cex.lab=1.4,cex.axis=1.3)
abline(v=2,col="red")

sum(t.samples>=2)/S


