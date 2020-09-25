# Tutorial 1 R code #

##############
# Question 1 #
##############

# 1(b)
theta <- seq(from=0,to=1,by=0.1)
y <- 57
n <- 100
likelihood <- dbinom(y, n, theta)
round(likelihood,4)

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(theta, likelihood, type="h", lwd=2.5,
    xlab=expression(theta), ylab=expression(paste("P(Y=57|",theta,")")))

# 1(c)
posterior <- likelihood/sum(likelihood)

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(theta, posterior,type="h",lwd=2.5,
    xlab=expression(theta), ylab=expression(paste("p(", theta, "|Y=57)")))


# 1(d)
theta2 <- seq(from=0, to=1, length.out=100)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(theta2, dbinom(y,n,theta2), type="l",lwd=2.5,
    xlab=expression(theta), ylab=expression(paste("p(", theta, ") p(Y=57| ", theta, ")")))


# 1(e)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(theta2, dbeta(theta2,58,44),type="l",lwd=2.5,xlab=expression(theta), 
    ylab=expression(paste("p(",theta,"|Y=57)")))



##############
# Question 2 #
##############

# 2(a)
a0 <- 11/12
b0 <- 11/6
theta <- seq(from=0,to=1, by=0.001)

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(theta, dbeta(theta,a0,b0), type="l",lwd=2.5,
    xlab=expression(theta), ylab=expression(paste("p(", theta, ")")))

beta_var <- function(a,b) {a*b/((a+b)^2*(a+b+1))}
beta_var(a0,b0)
sqrt(beta_var(a0,b0))


2(b)
y <- 3
n <- 11
(a <- a0+y)
(b <- b0+n-y)

par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(theta, dbeta(theta,a,b),type="l",lwd=2.5,
    xlab=expression(theta), ylab="Density")
points(theta, dbeta(theta,a0,b0),type="l",lwd=2.5,lty=2)
legend("topright",legend=c("prior","posterior"),lty=c(1,2), 
    lwd=2.5,cex=1)

(mean <- a/(a+b))
(median <- qbeta(0.5, a, b))
(mode <- (a-1)/(a+b-2))


2(c)
# HPD region
require(TeachingDemos)
(hpdr <- hpd(qbeta,shape1=a,shape2=b,conf=0.5))
hpdr[2] - hpdr[1]

# quantile-based CI
(CI <- qbeta(c(0.25,0.75),a,b))
CI[2] -CI[1]


2(d)
86/433



##############
# Question 3 #
##############

# 3(a)
a0 <- 1; b0 <- 1
y1 <- 2; n1 <- 15
(a <- a0+y1)
(b <- b0+n1-y1)


# 3(c)
n2 <- 278
y2 <- seq(from=0, to=n2, by=1)
predictive_density <- function(y2, n2, a, b){ 
 choose(n2,y2)*beta(y2+a,n2-y2+b)/beta(a,b) 
}
y2prob <- predictive_density(y2, n2, a, b)


par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(y2,y2prob,type="h",lwd=2.5,ylab="Density",
    xlab=expression(y[2]),main="predictive density")

# mean
(post.mean <- sum(y2*y2prob))
# standard deviation
post.var <- sum(y2^2*y2prob) - post.mean^2
(post.sd <- sqrt(post.var))


# 3(d)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.1,0.8,0))
plot(y2,dbinom(y2,n2,2/15),type="h",lwd=2.5,
    xlab=expression(y[2]),ylab="Density",main="Binomial(278,2/15)")

278*2/15
sqrt(278*2/15*13/15)




