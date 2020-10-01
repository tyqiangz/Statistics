###### Chapter 6 Rejection Sampling and Importance Sampling	######

### Example: Stomach cancer data

require(LearnBayes)            # install and load the R package "LearnBayes"
require(mvtnorm)			   # need "mvtnorm" for multivariate t density

data(cancermortality)          # load the dataset "cancermortality" 
y <- cancermortality$y
n <- cancermortality$n
NC <- length(y)                # no. of cities            

logpost <- function(theta, y, n){
	eta <- 1/(1+exp(-theta[1]))
	K <- exp(theta[2])
	(sum(lbeta(K*eta+y,K*(1-eta)+n-y) - lbeta(K*eta,K*(1-eta)))
     + theta[2] - 2*log(1+exp(theta[2])))
}

(out <- optim(par=c(-7,7.5),fn=logpost,hessian=TRUE,control=list(fnscale=-1),y=y,n=n))
(post.mode <- out$par)
(post.cov <- -solve(out$hessian))




## rejection sampling
data <- list(y=y,n=n,post.mode=post.mode,post.cov=post.cov,df=4)
diff <- function(theta, data){
    y <- data$y
    n <- data$n
    post.mode <- data$post.mode
    post.cov <- data$post.cov
    df <- data$df
    return(logpost(theta,y,n) - 
       dmvt(theta,delta=post.mode,sigma=2*post.cov,df=df))
}

# draw the contour plot on page 16
N <- 100
theta1 <- seq(from=-8, to=-4.5, length.out=N)
theta2 <- seq(from=3, to=17, length.out=N)
z <- matrix(0,N,N)
for (i in 1:N){
 for (j in 1:N){
  theta <- c(theta1[i], theta2[j])
  z[i,j] <- diff(theta, data)    
}}

par(mar=c(4,4.5,1,1))
contour(x=theta1,y=theta2,z,col="blue",nlevels=40,
    xlab=expression(theta[1]),ylab=expression(theta[2]),
    cex.axis=1.3,cex.lab=1.6,labcex=1.1)


# optimize the diff() function
(out1 <- optim(par=c(-6.8,12),fn=diff,
               control=list(fnscale=-1),data=data))

logM <- out1$value
data$logM <- logM

# rejection sampling main function
rej_sampling <- function(S, data){
    post.mode <- data$post.mode     
    df <- data$df
    post.cov <- data$post.cov        
    logM <- data$logM
    theta.samples <- matrix(0,S,length(post.mode))
    s <- 0             # no. of samples collected
    T <- 0             # no. of iterations performed
    while (s < S){
        T <- T+1
        theta <- rmvt(1,delta=post.mode,sigma=2*post.cov,df=df)
        U <- runif(1)
        if (log(U) < diff(theta,data) - logM){
            s <- s+1
            theta.samples[s,] <- theta
        } 
    }
    accept.rate <- S/T
    list(accept.rate=accept.rate, theta.samples=theta.samples)
}

# perform rejection sampling until S samples are collected
set.seed(1)
S <- 10000
output <- rej_sampling(S, data)
output$accept.rate

# draw the plot on page 20
N <- 100
theta1 <- seq(from=-8, to=-4.5, length.out=N)
theta2 <- seq(from=3, to=17, length.out=N)
z <- matrix(0,N,N)
for (i in 1:N){
 for (j in 1:N){
   theta <- c(theta1[i], theta2[j])
   z[i,j] <- logpost(theta,y,n)    
}}

par(mar=c(4,4,1,1))
par(mgp=c(2.3,1,0))
contour(x=theta1,y=theta2,z,col="blue",nlevels=50,
    xlab=expression(theta[1]),ylab=expression(theta[2]),
    cex.axis=1.4,cex.lab=1.8,drawlabels=FALSE)
points(output$theta.samples[,1],output$theta.samples[,2], col=rgb(1,0,0,alpha=0.2),pch=19)
      


## important sampling

logf <- function(theta2,theta1,y,n){
	J <- length(y)   # no. of cities
	eta <- 1/(1+exp(-theta1))
	K <- exp(theta2)
	L <- rep(0, length(theta2))
	for (j in 1:J){
		L <- L + lbeta(K*eta+y[j],K*(1-eta)+n[j]-y[j]) - lbeta(K*eta,K*(1-eta)) 
	}
	L <- L + theta2 - 2*log(1+K)
	return(L)
}


df <- 4
theta1 <- post.mode[1]
theta2.grid <- seq(from=1,to=20,by=0.1)

logw1 <- logf(theta2.grid,theta1,y,n) - dnorm(theta2.grid, mean=8,sd=2,log=TRUE)
w1 <- exp(logw1)

logw2 <- logf(theta2.grid,theta1,y,n) - dmvt(matrix(theta2.grid,ncol=1),delta=8,sigma=matrix(2,1,1),df=df)
w2 <- exp(logw2)


# density plots on page 28
par(mar=c(3.5,4.2,1.5,0.5))
par(mgp=c(2.5,0.8,0))
par(mfrow=c(1,4))
plot(theta2.grid, exp(logf(theta2.grid,theta1,y,n)), type="l",
  xlab=expression(theta[2]),ylab=expression(paste("f(",theta[2],")")),
  cex.lab=1.8,cex.axis=1.3)
plot(theta2.grid,dnorm(theta2.grid, mean=8, sd=2),type="l",ylim=c(0,0.27),
  xlab=expression(theta[2]),ylab="density",col="blue",cex.lab=1.8,cex.axis=1.3)
points(theta2.grid,
  dmvt(matrix(theta2.grid,ncol=1),delta=8,sigma=matrix(2,1,1),df=df,log=FALSE),
  type="l",col="red")
legend("topright", legend=c("N(8,4)", expression(paste(t[4],"(8,2)"))),
  lty=1,col=c("blue","red"),cex=1.2)
plot(theta2.grid,w1,type="l",xlab=expression(theta[2]),col="blue",cex.lab=1.8,cex.axis=1.3,
  ylab="unnormalized weights",main="N(8,4)",cex.main=1.6)
plot(theta2.grid,w2,type="l",xlab=expression(theta[2]),col="red",cex.lab=1.8,cex.axis=1.3,
   ylab="unnormalized weights",main=expression(paste(t[4],"(8,2)")),cex.main=1.6)


# estimate posterior mean of theta2 conditional on theta1
set.seed(1)
S <- 10000
df <- 4
theta1 <- post.mode[1]
theta2.samples <- rmvt(S,delta=8,sigma=matrix(2,1,1),df=4)
logw <- logf(theta2.samples,theta1,y,n) - dmvt(theta2.samples,delta=8,sigma=matrix(2,1,1),df=df)
w <- exp(logw - max(logw))
W <- w/sum(w)
(est <- sum(W*theta2.samples))
(se <- sqrt(sum((W*(theta2.samples-est))^2)))
par(mfrow=c(1,1))
plot(W)		

# estimate the marginal mean of theta2
# perform importance sampling for the joint posterior of (theta1,theta2)
set.seed(1)
S <- 10000
df <- 4
theta.samples <- rmvt(S, delta=post.mode, sigma=2*post.cov, df=df)
lf <- rep(0,S)
for (s in 1:S){ lf[s] <- logpost(theta.samples[s,],y,n) }
logw <- lf - dmvt(theta.samples,delta=post.mode,sigma=2*post.cov,df=df)
w <- exp(logw - max(logw))
W <- w/sum(w)
(est <- sum(W*theta.samples[,2]))
(se <- sqrt(sum((W*(theta.samples[,2] -est))^2)))

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
hist(W)		# the histogram on page 33



## Sampling importance resampling (SIR)
# a simple demonstration
set.seed(123)
N <- 10000
indices <- sample(1:S, size=N, prob=W, replace=TRUE)
theta.SIR <- theta.samples[indices,]	# theta.samples are generated from t above
# a plot of SIR samples
M <- 100
theta1 <- seq(from=-8, to=-4.5, length.out=M)
theta2 <- seq(from=3, to=17, length.out=M)
z <- matrix(0,M,M)
for (i in 1:M){
 for (j in 1:M){
   theta <- c(theta1[i], theta2[j])
   z[i,j] <- logpost(theta,y,n)    
}}
par(mar=c(4,4,1,1))
par(mgp=c(2.3,1,0))
contour(x=theta1,y=theta2,z,col="blue",nlevels=50,
    xlab=expression(theta[1]),ylab=expression(theta[2]),
    cex.axis=1.4,cex.lab=1.8,drawlabels=FALSE)
points(theta.SIR[,1],theta.SIR[,2], col=rgb(0.6,1,0.4,alpha=0.2),pch=19)
   


# leave out observation i
theta.samples <- theta.SIR
S <- nrow(theta.samples)
summary.full1 <- quantile(theta.samples[,1], c(0.05,0.5,0.95))
summary.full2 <- quantile(theta.samples[,2], c(0.05,0.5,0.95))
summary.obs1 <- matrix(0, 20, 3)
summary.obs2 <- matrix(0, 20, 3)
for (i in 1:20){
	eta <- 1/(1+exp(-theta.samples[,1]))
	K <- exp(theta.samples[,2])
	logw <- lbeta(K*eta,K*(1-eta)) - lbeta(K*eta+y[i],K*(1-eta)+n[i]-y[i])
	w <- exp(logw - max(logw))
	W <- w/sum(w)
	indices <- sample(1:S, size=N, prob=W, replace=TRUE)
	theta.SIR <- theta.samples[indices,]
	summary.obs1[i,] <- quantile(theta.SIR[,1], c(0.05,0.5,0.95))
	summary.obs2[i,] <- quantile(theta.SIR[,2], c(0.05,0.5,0.95))
}


par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
par(mfrow=c(1,2))
plot(c(0,0,0),summary.full1, xlim=c(-1,21),ylim=c(-7.4,-6.1),
  type="b",col="blue",lwd=2,xlab="Observation removed", 
  ylab=expression(theta[1]),cex.lab=1.2)
for (i in 1:20){
  lines(c(i,i,i),summary.obs1[i,],type="b")
}
lines(rep(15,3),summary.obs1[15,],type="b",col="red",lwd=2)
lines(rep(10,3),summary.obs1[10,],type="b",col="green",lwd=2)
lines(rep(19,3),summary.obs1[19,],type="b",col="green",lwd=2)
plot(c(0,0,0),summary.full2, xlim=c(-1,21),ylim=c(5.3,11.1),
  type="b",col="blue",lwd=2,xlab="Observation removed", 
  ylab=expression(theta[2]),cex.lab=1.2)
for (i in 1:20){
  lines(c(i,i,i),summary.obs2[i,],type="b")
}
lines(rep(15,3),summary.obs2[15,],type="b",col="red",lwd=2)
lines(rep(10,3),summary.obs2[10,],type="b",col="green",lwd=2)
lines(rep(19,3),summary.obs2[19,],type="b",col="green",lwd=2)

# show the three "influential" observations
c(y[15],n[15])
c(y[10],n[10])
c(y[19],n[19])


