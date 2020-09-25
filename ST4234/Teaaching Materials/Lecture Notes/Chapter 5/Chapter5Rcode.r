###### Chapter 5 Normal Approximation and Laplace Approximation	######


### Example: Genetic linkage model

# plot the exact posterior
Y <- c(125, 18, 20, 34)
post0 <- function(theta, Y) {
	(2+theta)^(Y[1])*(1-theta)^(Y[2]+Y[3])*theta^(Y[4])
}
(Const <- integrate(post0, 0, 1, Y=Y)$value)

pdf("truepost.pdf",width=4,height=3.5)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
theta.grid <- seq(from=0, to=1, by=0.001)
plot(theta.grid, post0(theta.grid,Y=Y)/Const, type="l", lwd=2, 
	 ylab=expression(paste("p(",theta,"|y)")), xlab=expression(theta),
	 main="True posterior density")
dev.off()  

# find the mode of p(theta|y) and the Hessian
(out <- optimize(post0, interval=c(0,1), Y=Y, maximum=TRUE))
theta.hat <- out$maximum
(var.theta <- 1/(Y[1]/(2+theta.hat)^2+(Y[2]+Y[3])/(1-theta.hat)^2+Y[4]/theta.hat^2))

# plot the normal approximation
pdf("normalapprox1.pdf",width=5,height=4)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
theta.grid <- seq(from=0, to=1, by=0.001)
plot(theta.grid, post0(theta.grid,Y=Y)/Const, type="l", lwd=2, 
	 ylab=expression(paste("p(",theta,"|y)")), xlab=expression(theta),
	 main="Normal approximation")
points(theta.grid, dnorm(theta.grid,mean=theta.hat,sd=sqrt(var.theta)),
		type="l", lwd=2, lty=2, col="red")
legend("topleft",legend = c(expression(paste("posterior p(", theta, "|y)")),
		"normal approximation"),
		col=c("black","red"), lty=c(1,2), lwd=2, cex=0.8)
dev.off() 

# repeat the process for Y=(14,0,1,5)
Y <- c(14, 0, 1, 5)
Const2 <- integrate(post0, 0, 1, Y=Y)$value
(out <- optimize(post0, interval=c(0,1), Y=Y, maximum=TRUE))
theta.hat <- out$maximum
(var.theta <- 1/(Y[1]/(2+theta.hat)^2+(Y[2]+Y[3])/(1-theta.hat)^2+Y[4]/theta.hat^2))

pdf("normalapprox2.pdf",width=5,height=4)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
theta.grid <- seq(from=0, to=1, by=0.001)
plot(theta.grid, post0(theta.grid,Y=Y)/Const2, type="l", lwd=2, 
	 ylab=expression(paste("p(",theta,"|y)")), xlab=expression(theta),
	 main="Normal approximation")
points(theta.grid, dnorm(theta.grid,mean=theta.hat,sd=sqrt(var.theta)),
		type="l", lwd=2, lty=2, col="red")
legend("topleft",legend = c(expression(paste("posterior p(", theta, "|y)")),
		"normal approximation"),
		col=c("black","red"), lty=c(1,2), lwd=2, cex=0.8)
dev.off() 

## Laplace method
# calculate the Laplace approximation to denominator
Y <- c(125, 18, 20, 34)
# Y <- c(14, 0, 1, 5)
post0 <- function(theta, Y) {
	(2+theta)^(Y[1])*(1-theta)^(Y[2]+Y[3])*theta^(Y[4])
}
Const <- integrate(post0, 0, 1, Y=Y)$value
out <- optimize(post0, interval=c(0,1), Y=Y, maximum=TRUE)
theta.hat <- out$maximum
s.hat <- 1/sqrt(Y[1]/(2+theta.hat)^2+
          (Y[2]+Y[3])/(1-theta.hat)^2+Y[4]/theta.hat^2)

# plot the Laplace approximation together with normal approximation
pdf("Laplace1.pdf",width=5,height=4)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
theta.grid <- seq(from=0, to=1, by=0.001)
plot(theta.grid, post0(theta.grid,Y=Y)/Const, type="l", lwd=2, 
	 ylab=expression(paste("p(",theta,"|y)")), xlab=expression(theta),
	 main="Laplace approximation")
points(theta.grid, dnorm(theta.grid,mean=theta.hat,sd=s.hat),
		type="l", lwd=2, lty=2, col="red")
points(theta.grid, post0(theta.grid,Y)/(sqrt(2*pi)*s.hat*post0(theta.hat,Y)),
		type="l", lwd=3, lty=3, col="blue")		
legend("topleft",legend = c(expression(paste("posterior p(", theta, "|y)")),
		"normal approximation", "Laplace approximation"),
		col=c("black","red","blue"), lty=c(1,2,3), lwd=2, cex=0.8)
dev.off() 		  

# estimate the posterior mean of theta
numerator <- function(theta, Y) {
	theta * (2+theta)^(Y[1])*(1-theta)^(Y[2]+Y[3])*theta^(Y[4])
}
(integrate(numerator, 0, 1, Y=Y)$value)/Const      # true posterior mean
(theta.hat)      # Approximation 1

out.star <- optimize(numerator, interval=c(0,1), Y=Y, maximum=TRUE)
theta.star <- out.star$maximum
s.star <- 1/sqrt(1/theta.star^2 + Y[1]/(2+theta.star)^2+
          (Y[2]+Y[3])/(1-theta.star)^2 + Y[4]/theta.star^2)
((s.star*theta.star*post0(theta.star,Y))/
 (s.hat*post0(theta.hat,Y)))      # Approximation 2



 
 

### Example: Stomach cancer data

# contour plot of cancer mortality dataset
require(LearnBayes)            # install and load the R package "LearnBayes"
data(cancermortality)          # load the dataset "cancermortality" 
y <- cancermortality$y
n <- cancermortality$n
NC <- length(y)                # no. of cities            

logpost0 <- function(eta,K,y,n){
	(sum(lbeta(K*eta+y,K*(1-eta)+n-y) - lbeta(K*eta,K*(1-eta)))
     -log(eta) - log(1-eta) -2*log(1+K))
}

N <- 100
K <- seq(from=1, to=20000, length.out=N)
eta <- seq(from=0, to=0.003, length.out=N)
z <- matrix(0,N,N)
for (i in 1:N){
	for (j in 1:N){
		z[i,j] <- logpost0(eta[i],K[j],y,n)    
}}

pdf("contour01.pdf",width=4,height=3.5)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
contour(x=eta,y=K,z,col="blue",nlevels=40, 
    xlab=expression(eta),ylab="K",cex.lab=1.5,cex.axis=1.2)
dev.off()      


# contour plot of transformed parameters

logpost <- function(theta,y,n){
	eta <- 1/(1+exp(-theta[1]))
	K <- exp(theta[2])
	(sum(lbeta(K*eta+y,K*(1-eta)+n-y) - lbeta(K*eta,K*(1-eta)))
     + theta[2] - 2*log(1+exp(theta[2])))
}

N <- 100
theta1 <- seq(from=-8, to=-4.5, length.out=N)
theta2 <- seq(from=3, to=17, length.out=N)
z <- matrix(0,N,N)
for (i in 1:N){
	for (j in 1:N){
		theta <- c(theta1[i], theta2[j])
		z[i,j] <- logpost(theta,y,n)    
}}


pdf("contour1.pdf",width=4,height=3.5)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
contour(x=theta1,y=theta2,z,col="blue",nlevels=80,
    xlab=expression(theta[1]),ylab=expression(theta[2]),
    cex.axis=1.2,cex.lab=1.5,drawlabels=FALSE)
dev.off()      



# normal approximation
(out <- optim(par=c(-7,7.5),fn=logpost,hessian=TRUE,control=list(fnscale=-1),y=y,n=n))
(post.mode <- out$par)
(post.cov <- -solve(out$hessian))

require(mvtnorm)
N <- 100
theta1 <- seq(from=-8, to=-4.5, length.out=N)
theta2 <- seq(from=3, to=17, length.out=N)
z <- matrix(0,N,N)
for (i in 1:N){
	for (j in 1:N){
		theta <- c(theta1[i], theta2[j])
		z[i,j] <- dmvnorm(theta,mean=post.mode,sigma=post.cov)    
}}

pdf("contour2.pdf",width=4,height=3.5)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
contour(x=theta1,y=theta2,z,col="blue",drawlabels=FALSE,
    xlab=expression(theta[1]),ylab=expression(theta[2]),
    cex.axis=1.1,cex.lab=1.3)
dev.off()


post.sd <- sqrt(diag(post.cov))
qnorm(c(0.05,0.95),mean=post.mode[1],sd=post.sd[1])
qnorm(c(0.05,0.95),mean=post.mode[2],sd=post.sd[2])


# marginal approximation plot
# for eta
plot(density())






